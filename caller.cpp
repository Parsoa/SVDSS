#include "caller.hpp"

void Caller::run() {
    config = Configuration::getInstance();

    load_chromosomes(config->reference);
    cout << "Loaded all chromosomes" << endl ;
    load_input_sfs() ;

    ovcf.open(config->workdir + "/svs_hybrid.vcf");
    osam.open(config->workdir + "/poa.sam");
    // --- SAM header
    osam << "@HD\tVN:1.4" << endl;
    for (int i = 0; i < chromosomes.size(); ++i) {
        osam << "@SQ\tSN:" << chromosomes[i] << "\t" << "LN:" << strlen(chromosome_seqs[chromosomes[i]]) << endl ;
    }
    // --- VCF header
    print_vcf_header() ;

    //vector<vector<Consensus>> alignments(config->threads) ; // produce a SAM file of consensus alignments
    //vector<vector<SV>> svs(config->threads) ; // produce a SAM file of consensus alignments
    //#pragma omp parallel for num_threads(config->threads) schedule(static,1)
    //for(int i = 0; i < chromosomes.size(); i++) {
    //    string chrom = chromosomes[i] ;
    //    int t = i % config->threads ;
    //    cout << "Processing chromosome " << chrom << ".. " << endl ;

    //    Extender extender = Extender(chrom, &SFSs) ;
    //    extender.run(4) ;
    //    alignments[t].insert(alignments[t].begin(), extender.alignments.begin(), extender.alignments.end()) ;
    //    svs[t].insert(svs[t].begin(), extender.svs.begin(), extender.svs.end()) ;

    //    cout << svs[t].size() << " SVs." << endl ;

    //    //Clipper clipper(chrom, extender.clips);
    //    //clipper.call(reference[chrom], sv_tree);
    //    //svs[omp_get_thread_num()].insert(svs[omp_get_thread_num()].begin(), clipper.svs.begin(), clipper.svs.end()) ;
    //}
    //for (int i = 0; i < config->threads; i++) {
    //    for (int j = 0; j < alignments[i].size(); j++) {
    //        const auto& c = alignments[i][j] ;
    //        osam << c.chrom << ":" << c.s + 1 << "-" << c.e + 1 << "\t"
    //            << "0"
    //            << "\t" << c.chrom << "\t" << c.s + 1 << "\t"
    //            << "60"
    //            << "\t" << c.cigar << "\t"
    //            << "*"
    //            << "\t"
    //            << "0"
    //            << "\t"
    //            << "0"
    //            << "\t" << c.seq << "\t"
    //            << "*" << endl ;
    //    }
    //    for (const SV& sv: svs[i]) {
    //        ovcf << sv << endl ;
    //    }
    //}


    Extender extender = Extender(&SFSs) ;
    extender.run(config->threads) ;
    //Clipper clipper(chrom, extender.clips);
    //clipper.call(reference[chrom], sv_tree);
    //svs[omp_get_thread_num()].insert(svs[omp_get_thread_num()].begin(), clipper.svs.begin(), clipper.svs.end()) ;

    for (int j = 0; j < extender.alignments.size(); j++) {
        const auto& c = extender.alignments[j] ;
        osam << c.chrom << ":" << c.s + 1 << "-" << c.e + 1 << "\t"
            << "0"
            << "\t" << c.chrom << "\t" << c.s + 1 << "\t"
            << "60"
            << "\t" << c.cigar << "\t"
            << "*"
            << "\t"
            << "0"
            << "\t"
            << "0"
            << "\t" << c.seq << "\t"
            << "*" << endl ;
    }
    
    std::sort(extender.svs.begin(), extender.svs.end()) ;
    for (const SV& sv: extender.svs) {
        ovcf << sv << endl ;
    }

    osam.close() ;
    ovcf.close() ;
}

void Caller::load_input_sfs() {
    int threads = config->threads ; 
    int num_batches = config->aggregate_batches ;
    int num_threads = num_batches < threads ? num_batches : threads ;
    vector<unordered_map<string, vector<SFS>>> _SFSs(num_batches) ;
    cout << "Loading assmbled SFS.." << endl ;
    #pragma omp parallel for num_threads(num_threads)
    for (int j = 0; j < num_batches; j++) {
        string s_j = std::to_string(j) ;
        string inpath = config->workdir + "/solution_batch_" + s_j + ".assembled.sfs" ;
        cout << "[I] Loading SFS from " << inpath << endl ;
        ifstream inf(inpath) ;
        string line ;
        if (inf.is_open()) {
            string info[4];
            string read_name;
            while (getline(inf, line)) {
                stringstream ssin(line);
                int i = 0;
                while (ssin.good() && i < 4) {
                    ssin >> info[i++];
                }
                if (info[0].compare("*") != 0) {
                    read_name = info[0];
                    _SFSs[j][read_name] = vector<SFS>();
                }
                _SFSs[j][read_name].push_back(SFS(stoi(info[1]), stoi(info[2]), stoi(info[3]), true)) ;
            }
        }
    }
    int r = 0 ;
    int c = 0 ;
    for (int j = 0; j < num_batches; j++) {
        lprint({"Batch", to_string(j), "with", to_string(_SFSs[j].size()), "strings."});
        r += _SFSs[j].size() ;
        SFSs.insert(_SFSs[j].begin(), _SFSs[j].end()) ;
        for (auto& read: _SFSs[j]) {
            c += read.second.size() ;
        }
    }
    lprint({"Loaded", to_string(c), "SFS strings on", to_string(r), " reads."}) ;
}

void Caller::print_vcf_header() {
    ovcf << "##fileformat=VCFv4.2" << endl;
    ovcf << "##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz" << endl;
    for (int i = 0; i < chromosomes.size(); ++i) {
        ovcf << "##contig=<ID=" << chromosomes[i] << ",length=" << strlen(chromosome_seqs[chromosomes[i]]) << ">" << endl;
    }
    ovcf << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << endl;
    ovcf << "##INFO=<ID=VARTYPE,Number=A,Type=String,Description=\"Variant class\">" << endl;
    ovcf << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant type\">" << endl;
    ovcf << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
    ovcf << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
    ovcf << "##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description=\"Number of alignments supporting this record\">" << endl;
    ovcf << "##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Total number of alignments covering this locus\">" << endl;
    ovcf << "##INFO=<ID=ASCORE,Number=1,Type=Integer,Description=\"Alignment score\">" << endl;
    ovcf << "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of variations on same consensus\">" << endl;
    ovcf << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << endl;
    ovcf << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    ovcf << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEFAULT" << endl;
}
