#include <sys/stat.h>
#include <sys/types.h>

#include "htslib/sam.h"
#include "htslib/hts.h"

#include "haplotyper.hpp"

Haplotyper::Haplotyper() {
    config = Configuration::getInstance() ;
}

vector<SV> Haplotyper::assemble_reads(Cluster cluster) {
    struct stat info ;
    string name = cluster.chrom + "_" + to_string(cluster.s) + "_" + to_string(cluster.e) ;
    string cluster_id = cluster.get_id() ;
    // make working directory
    string dir = config->workdir + "/miniasm" ;
    if (stat(dir.c_str(), &info) != 0) {
        int check = mkdir(dir.c_str(), 0777) ;
        if (check != 0) {
            cerr << "Error creating output directory " << dir << ".." << endl ;
            exit(check) ;
        }
    }
    dir += "/" + cluster_id ;
    if (stat(dir.c_str(), &info) != 0) {
        int check = mkdir(dir.c_str(), 0777) ;
        if (check != 0) {
            cerr << "Error creating output directory " << dir << ".." << endl ;
            exit(check) ;
        }
    }
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
    lprint({"Generated", std::to_string(unitigs.size()), "unitigs."}) ;
    vector<SV> svs ;
    if (unitigs.size() == 0) {
        lprint({"Error assembling cluster", cluster_id}, 'E') ;
    } else {
        // create a FASTA file for unitigs
        string unitig_fasta = dir + "/unitigs.fasta" ;
        //cout << unitig_fasta << endl ;
        ofstream unitig_fasta_file(unitig_fasta) ;
        unitig_fasta_file << ">" + name << endl ;
        unitig_fasta_file << unitigs[0] << endl ;
        unitig_fasta_file.close() ; 
        // get reference sequence corresponding to cluster
        int s = cluster.s ;
        int e = cluster.e ;
        string ref_fasta = dir + "/ref.fasta" ;
        int d = 2 * unitigs[0].length() ;
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
        samFile *bam_file = hts_open(bam_path.c_str(), "r") ;
        if (bam_file == nullptr) {
            lprint({"Error openning consesnsus alignemnt for", name}, 'E') ;
            return svs ;
        }
        bam_hdr_t *bam_header = sam_hdr_read(bam_file) ;
        bam1_t *aln = bam_init1();
        cout << "Generated " << unitigs.size() << " unitigs.." << endl ;
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
                    string alt_allele = ref_allele + unitigs[0].substr(pos, cigar[i].first) ;
                    sv = SV("INS", cluster.chrom, s - d + ref_pos - 1, ref_allele, alt_allele, cigar[i].first, 1, -1, -1) ;
                    sv.idx = sv.idx + "_" + cluster.get_id() ;
                    pos += cigar[i].first ;
                }
                if (cigar[i].second == BAM_CDEL) {
                    string ref_allele = ref_seq.substr(ref_pos - 1, cigar[i].first + 1) ;
                    string alt_allele = ref_allele.substr(0, 1) ;
                    sv = SV("DEL", cluster.chrom, s - d + ref_pos - 1, ref_allele, alt_allele, cigar[i].first, 1, -1, -1) ;
                    sv.idx = sv.idx + "_" + cluster.get_id() ;
                    ref_pos += cigar[i].first ;
                }
                if (cigar[i].second == BAM_CMATCH) {
                    pos += cigar[i].first ;
                    ref_pos += cigar[i].first ;
                }
                if (abs(sv.l) > 10) {
                    svs.push_back(sv) ;
                }
            }
        }
        sam_close(bam_file) ;
    }
    for (auto sv: svs) {
       cout << "[S]" << sv << endl ; 
    }
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
