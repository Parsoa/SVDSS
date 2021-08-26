#include "caller.hpp"

void Caller::run() {
    config = Configuration::getInstance();

    load_chromosomes(config->reference);
    cout << "Loaded all chromosomes" << endl ;

    ovcf.open(config->workdir + "/svs.vcf");
    osam.open(config->workdir + "/poa.sam");

    map<string, vector<SV>> clips ; // this will contain variations from clipping
    map<string, vector<SV>> indels ; // this will contain all variations from Is and Ds

    read_bam = hts_open(config->bam.c_str(), "r") ;
    read_bamhdr = sam_hdr_read(read_bam) ;
    read_bamindex = sam_index_load(read_bam, config->bam.c_str()) ;
    osam << read_bamhdr->text ;
    
    for (uint i = 0; i < chromosomes.size(); ++i) {
        string chrom = chromosomes[i];
        cout << "Processing chromosome " << chrom << ".. " << endl ;

        Insdeller indeler = Insdeller(chrom) ;
        indeler.call(chromosome_seqs[chrom], osam) ;
        indels[chrom] = indeler.osvs ;

        //Clipper clipper = Clipper(chrom) ;
        //clipper.call(chromosome_seqs[chrom], indeler.vartree);
        //clipper.call("", indeler.vartree);
        //clips[chrom] = clipper.indels;
    }

    print_vcf_header(chromosome_seqs, ovcf);
    for(const string &chrom : chromosomes) {
        for (const SV &sv : indels[chrom]) {
            ovcf << sv << endl;
        }
        for (const SV &sv : clips[chrom]) { 
            ovcf << sv << endl;
        }
    }

    osam.close();
    ovcf.close();
}
