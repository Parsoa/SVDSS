#include "caller.hpp"

void Caller::run()
{
    config = Configuration::getInstance();

    load_chromosomes(config->reference);

    ovcf.open(config->workdir + "/svs.vcf");
    osam.open(config->workdir + "/poa.sam");

    sfs_bam = hts_open(config->sfsbam.c_str(), "r");
    sfs_bamhdr = sam_hdr_read(sfs_bam);
    sfs_bamindex = sam_index_load(sfs_bam, config->sfsbam.c_str());

    read_bam = hts_open(config->bam.c_str(), "r");
    read_bamhdr = sam_hdr_read(read_bam);
    read_bamindex = sam_index_load(read_bam, config->bam.c_str());

    map<string, list<SV>> osvs; // this will contain all variations from Is and Ds
    map<string, list<SV>> oimprsvs; // this will contain variations from clipping
    // FIXME: we can easily parallelize this for but dumping to osam will break I think
    osam << read_bamhdr->text;
    for (uint i = 0; i < chromosomes.size(); ++i)
    {
        string chrom = chromosomes[i];

        Insdeller idler = Insdeller(chrom, sfs_bam, sfs_bamhdr, sfs_bamindex, read_bam, read_bamhdr, read_bamindex);
        idler.call(chromosome_seqs[chrom], osam);
        osvs[chrom] = idler.osvs;

        Clipler cler = Clipler(chrom, sfs_bam, sfs_bamhdr, sfs_bamindex);
        cler.call(chromosome_seqs[chrom], idler.vartree);
        oimprsvs[chrom] = cler.osvs;
    }

    print_vcf_header(chromosome_seqs, ovcf);
    for(const string &chrom : chromosomes)
    {
        for (const SV &sv : osvs[chrom]) 
            ovcf << sv << endl;
        for (const SV &sv : oimprsvs[chrom]) 
            ovcf << sv << endl;
    }

    osam.close();
    ovcf.close();
}