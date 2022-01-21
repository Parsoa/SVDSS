#include "vcf.hpp"

unordered_map<string, vector<vcf_variant_t>> load_vcf_file(string path) {
    std::unordered_map<std::string, std::vector<vcf_variant_t>> vcf_variants ;
    std::ifstream vcf_file(path) ;
    std::string line ;
    unordered_map<string, int> header ;
    int i = 0 ;
    int n = 0 ;
    string chrom = "UNKNOWN" ;
    while (std::getline(vcf_file, line)) {
        if (line[0] == '#' and line[1] == '#') {
            continue ;
        } else if (line[0] == '#') {
            //TODO parse header
            istringstream iss(line) ;
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
        } else {
            istringstream iss(line) ;
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
            if (chrom != tokens[0]) {
                chrom = tokens[0] ;
            }
            vcf_variant_t vcf_variant ;
            vcf_variant.chrom = tokens[0] ;
            vcf_variant.pos = std::stoi(tokens[1]) ;
            vcf_variant.ref = tokens[3] ;
            int l = vcf_variants[chrom].size() ;
            if (l != 0 && vcf_variants[chrom][l - 1] == vcf_variant) {
                vcf_variants[chrom][l - 1].alleles[1] = tokens[4] ;
                vcf_variants[chrom][l - 1].svlen = max(vcf_variants[chrom][l - 1].svlen, int(tokens[4].length() - tokens[3].length())) ;
            } else {
                vcf_variant.alleles[0] = tokens[4] ;
                vcf_variant.alleles[1] = "$" ;
                vcf_variants[chrom].push_back(vcf_variant) ;
                vcf_variants[chrom][l - 1].svlen = tokens[4].length() - tokens[3].length() ;
            }
            n++ ;
        }
    }
    lprint({"Loaded", to_string(n), "variants from", path});
    return vcf_variants ;
}


void print_vcf_header(const unordered_map<string,  char*> &ref, ofstream &o, const string &sample)
{
    // TODO: improve and fix
    o << "##fileformat=VCFv4.2" << endl;
    // print("##fileDate=", date.today().strftime("%Y%m%d"), sep="")
    o << "##reference=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/hg38.no_alt.fa.gz" << endl;
    for (unordered_map<string,  char*>::const_iterator it = ref.begin(); it != ref.end(); ++it)
        o << "##contig=<ID=" << it->first << ",length=" << strlen(it->second) << ">" << endl;
    o << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << endl;
    o << "##INFO=<ID=VARTYPE,Number=A,Type=String,Description=\"Variant class\">" << endl;
    o << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Variant type\">" << endl;
    o << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
    o << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << endl;
    o << "##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description=\"Number of alignments supporting this record\">" << endl;
    o << "##INFO=<ID=COV,Number=1,Type=Integer,Description=\"Total number of alignments covering this locus\">" << endl;
    o << "##INFO=<ID=AS,Number=1,Type=Integer,Description=\"Alignment score\">" << endl;
    o << "##INFO=<ID=NV,Number=1,Type=Integer,Description=\"Number of variations on same consensus\">" << endl;
    o << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">" << endl;
    o << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    o << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample << endl;
}