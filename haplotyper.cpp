#include <sys/stat.h>
#include <sys/types.h>

#include "htslib/sam.h"
#include "htslib/hts.h"

#include "haplotyper.hpp"

Haplotyper::Haplotyper() {
    config = Configuration::getInstance() ;
}

vector<SV> Haplotyper::map_to_chm13(Cluster cluster, string assembly) {
    string ref_path = "/share/hormozdiarilab/Codes/Stella/output/T2T/v1.1/chr21.fasta" ; 
    string ref_seq = load_chromosome(ref_path) ;
    cout << ref_seq.length() << endl ;
    string cluster_id = cluster.chrom + "_" + to_string(cluster.s) + "_" + to_string(cluster.e) ;
    int check = make_working_dir(cluster) ;
    if (check != 0) {
        exit(check) ;
    }
    string dir = config->workdir + "/miniasm/" + cluster_id ;
    string poa_fasta = dir + "/poa.fasta" ;
    cout << poa_fasta << endl ;
    ofstream poa_fasta_file(poa_fasta) ;
    poa_fasta_file << ">" + cluster_id << endl ;
    poa_fasta_file << assembly << endl ;
    poa_fasta_file.close() ; 
    // map with minimap
    string bam_path = dir + "/aln.bam" ;
    string command = "minimap2 -ax map-pb -t 1 " + ref_path + " " + poa_fasta + " 2>/dev/null | samtools view -b > " + bam_path ;
    cout << command << endl ;
    system(command.c_str()) ; 
    auto svs = call_svs(cluster, bam_path, ref_seq, assembly, 0) ;
    return svs ;
}

int Haplotyper::make_working_dir(Cluster cluster) {
    struct stat info ;
    string cluster_id = cluster.chrom + "_" + to_string(cluster.s) + "_" + to_string(cluster.e) ;
    // make working directory
    string dir = config->workdir + "/miniasm" ;
    if (stat(dir.c_str(), &info) != 0) {
        int check = mkdir(dir.c_str(), 0777) ;
        if (check != 0) {
            cerr << "Error creating output directory " << dir << ".." << endl ;
            return check ;
        }
    }
    dir += "/" + cluster_id ;
    if (stat(dir.c_str(), &info) != 0) {
        int check = mkdir(dir.c_str(), 0777) ;
        if (check != 0) {
            cerr << "Error creating output directory " << dir << ".." << endl ;
            return check ;
        }
    }
    return 0 ;
}

vector<SV> Haplotyper::assemble_reads(Cluster cluster) {
    string cluster_id = cluster.chrom + "_" + to_string(cluster.s) + "_" + to_string(cluster.e) ;
    int check = make_working_dir(cluster) ;
    if (check != 0) {
        exit(check) ;
    }
    string dir = config->workdir + "/miniasm/" + cluster_id ;
    auto fastq_path = dir + "/reads.fastq" ;
    ofstream fastq_file(fastq_path) ;
    for (const auto& fragment: cluster.fragments) {
        fastq_file << "@" + fragment.name << endl ;
        fastq_file << fragment.seq << endl ;
        fastq_file << "+" << endl ;
        fastq_file << fragment.qual << endl ;
    }
    fastq_file.close() ;
    // generate overlap graph
    string graph_path = dir + "/reads.paf.gz" ;
    string command = "minimap2 -x ava-pb -t 1 " + fastq_path + " " + fastq_path + " 2>/dev/null | gzip -1 > " + graph_path ; 
    //lprint({"Running", command}) ;
    system(command.c_str()) ;
    // run miniasm
    string assembly_path = dir + "/reads.gfa" ;
    command = "miniasm -e 0 -n 0 -c 1 -1 -2 -s 200 -f " + fastq_path + " " + graph_path +  " > " + assembly_path + " 2>/dev/null";
    //lprint({"Running", command}) ;
    system(command.c_str()) ;
    // load assembly file
    ifstream gfa_file(assembly_path) ;
    string line ;
    vector<string> unitigs ;
    while (getline(gfa_file, line)) {
        istringstream iss(line) ;
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
        unitigs.push_back(tokens[2]) ;
    }
    gfa_file.close() ;
    vector<SV> svs ;
    lprint({"Generated", std::to_string(unitigs.size()), "unitigs."}) ;
    if (unitigs.size() == 0) {
        lprint({"Error assembling cluster", cluster.get_id()}, 'E') ;
    } else {
        // 
        //auto svs_ = map_to_chm13(cluster, unitigs[0]) ;
        //if (svs_.size() != 0) {
        //    cout << "[E]" << " assembly may be invalid." << endl ;
        //    int a ; 
        //    cin >> a ;
        //}
        // create a FASTA file for unitigs
        string unitig_fasta = dir + "/unitigs.fasta" ;
        //cout << unitig_fasta << endl ;
        ofstream unitig_fasta_file(unitig_fasta) ;
        unitig_fasta_file << ">" + cluster_id << endl ;
        unitig_fasta_file << unitigs[0] << endl ;
        unitig_fasta_file.close() ; 
        // get reference sequence corresponding to cluster
        int s = cluster.s ;
        int e = cluster.e ;
        string ref_fasta = dir + "/ref.fasta" ;
        int d = 2 * unitigs[0].length() ;
        cout << "Generated " << unitigs.size() << " unitigs.." << endl ;
        //cout << s - d << " " << e - s + d << " " << strlen(chromosome_seqs[cluster.chrom]) << endl ; 
        string ref_seq = string(chromosome_seqs[cluster.chrom]).substr(s - d, (e - s) + 2 * d) ; 
        ofstream ref_fasta_file(ref_fasta) ;
        ref_fasta_file << ">" << cluster.chrom << endl ;
        ref_fasta_file << ref_seq << endl ;
        ref_fasta_file.close() ;
        // call SV from this, map to reference again
        string bam_path = dir + "/aln.bam" ;
        command = "minimap2 -ax map-pb -t 1 " + ref_fasta + " " + unitig_fasta + " 2>/dev/null | samtools view -b > " + bam_path ;
        //cout << command << endl ;
        system(command.c_str()) ; 
        // read BAM and CIGAR 
        svs = call_svs(cluster, bam_path, ref_seq, unitigs[0], s - d) ;
    }
    for (auto sv: svs) {
       cout << "[S]" << sv << endl ; 
    }
    return svs ;
}

vector<SV> Haplotyper::call_svs(Cluster cluster, string bam_path, string ref_seq, string query, int ref_start) {
    vector<SV> svs ;
    samFile *bam_file = hts_open(bam_path.c_str(), "r") ;
    if (bam_file == nullptr) {
        lprint({"Error openning consesnsus alignemnt for", cluster.get_id()}, 'E') ;
        return svs ;
    }
    bam_hdr_t *bam_header = sam_hdr_read(bam_file) ;
    bam1_t *aln = bam_init1();
    while (sam_read1(bam_file, bam_header, aln) >= 0) {
        if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY) {
            continue ;
        }
        // get CIGAR and process
        uint pos = 0 ;
        uint ref_pos = aln->core.pos ;
        auto cigar = decode_cigar(aln) ;
        for (uint i = 0; i < cigar.size(); i++) {
            SV sv ;
            if (cigar[i].second == BAM_CSOFT_CLIP) {
                pos += cigar[i].first;
            }
            if (cigar[i].second == BAM_CINS) {
                string ref_allele = ref_seq.substr(ref_pos - 1, 1) ;
                //cout << "INS " << pos << " " << unitigs[0].length() << " " << cigar[i].first << endl ;
                string alt_allele = ref_allele + query.substr(pos, cigar[i].first) ;
                sv = SV("INS", cluster.chrom, ref_start + ref_pos - 1, ref_allele, alt_allele, cigar[i].first, 1, -1, -1) ;
                sv.idx = sv.idx + "_" + cluster.get_id() ;
                pos += cigar[i].first ;
            }
            if (cigar[i].second == BAM_CDEL) {
                string ref_allele = ref_seq.substr(ref_pos - 1, cigar[i].first + 1) ;
                string alt_allele = ref_allele.substr(0, 1) ;
                sv = SV("DEL", cluster.chrom, ref_start + ref_pos - 1, ref_allele, alt_allele, cigar[i].first, 1, -1, -1) ;
                sv.idx = sv.idx + "_" + cluster.get_id() ;
                ref_pos += cigar[i].first ;
            }
            if (cigar[i].second == BAM_CMATCH) {
                pos += cigar[i].first ;
                ref_pos += cigar[i].first ;
            }
            if (abs(sv.l) >= 10) {
                svs.push_back(sv) ;
            }
        }
    }
    sam_close(bam_file) ;
    return svs ;
}

vector<SV> Haplotyper::haplotype(Cluster cluster) {
    if (cluster.fragments.size() > 200) {
        lprint({"Skipping cluster ", cluster.get_id(), "with", to_string(cluster.fragments.size()), "fragments."}, 'E') ;
        return {} ;
    }
    cout << "[C] Cluster " << cluster.get_id() << " with " << cluster.fragments.size() << " fragments.." << endl ;
    return assemble_reads(cluster) ;
}
