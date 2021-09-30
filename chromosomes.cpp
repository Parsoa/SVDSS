#include "chromosomes.hpp"

KSEQ_INIT(gzFile, gzread)

std::vector<std::string> chromosomes ;
std::unordered_map<std::string, char*> chromosome_seqs ;

string reverse_complement(string s) {
    string rc(s) ;
    int l = s.length() ;
    for (int i = 0; i < l; i++) {
        if (s[i] == 'A') {
            rc[l - 1 - i] = 'T' ;
        }
        if (s[i] == 'C') {
            rc[l - 1 - i] = 'G' ;
        }
        if (s[i] == 'G') {
            rc[l - 1 - i] = 'C' ;
        }
        if (s[i] == 'T') {
            rc[l - 1 - i] = 'A' ;
        }
    }
    return rc ;
}

string canonicalize(string s) {
    auto rc = reverse_complement(s) ;
    return rc < s ? rc : s ;
}

int get_reference_size(ifstream &fasta_file) {
    fasta_file.seekg(0, ios_base::end) ;
    int l = fasta_file.tellg() ;
    fasta_file.seekg(0, ios_base::beg) ;
    return l ;
}

void load_chromosomes(string path) {
    lprint({"Loading reference genome from", path + ".."});
    gzFile fp = gzopen(path.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        string name(seq->name.s) ;
        if (name.find('_') == -1) {
            lprint({"Extracted", seq->name.s, "with", to_string(l), "bases."});
            for (uint i = 0; i < l; i++) {
                seq->seq.s[i] = toupper(seq->seq.s[i]) ;
            }
            chromosomes.push_back(seq->name.s) ;
            char* s = (char*) malloc(sizeof(char) * (l + 1)) ;
            memcpy(s, seq->seq.s, l + 1) ;
            s[l] = '\0' ;
            chromosome_seqs[seq->name.s] = s;
        }
    }
    kseq_destroy(seq);
    gzclose(fp);
}

string load_chromosome(string path) {
    lprint({"Loading chromosome genome from", path, ".."}) ;
    gzFile fp = gzopen(path.c_str(), "r") ;
    kseq_t *seq = kseq_init(fp) ;
    int l ;
    string chrom ;
    while ((l = kseq_read(seq)) >= 0) {
        lprint({"Extracted", seq->name.s, "with", to_string(l), "bases"});
        for (uint i = 0; i < l; i++) {
            seq->seq.s[i] = toupper(seq->seq.s[i]) ;
        }
        char* s = (char*) malloc(sizeof(char) * (l + 1)) ;
        memcpy(s, seq->seq.s, l + 1) ;
        chrom = string(s) ;
        s[l] = '\0' ;
    }
    kseq_destroy(seq) ;
    gzclose(fp) ;
    return chrom ;
}
