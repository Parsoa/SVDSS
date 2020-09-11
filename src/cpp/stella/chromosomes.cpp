#include <iomanip>

#include "chromosomes.hpp"

using namespace std ;

std::vector<std::string> chromosomes ;
std::unordered_map<std::string, char*> chromosome_seqs ;

string canonicalize(string s) {
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
    return rc < s ? rc : s ;
}

int get_reference_size(ifstream &fasta_file) {
    fasta_file.seekg(0, ios_base::end) ;
    int l = fasta_file.tellg() ;
    fasta_file.seekg(0, ios_base::beg) ;
    return l ;
}

void load_chromosomes(string path) {
    cout << "Loading reference genome from " << path << ".." << endl ;
    ifstream fasta_file ;
    fasta_file.open(path, ios::binary) ;
    // maximum size of a chromosome, kinda arbitrary
    char* buffer = (char*) malloc(sizeof(char) * 300000000) ;
    int state ;
    uint64_t n = 0 ;
    std::string line ;
    std::getline(fasta_file, line) ;
    while (true) {
        if (line.substr(0, 4) == ">chr") {
            int l = line.length() ;
                if (line[4] == 'X' || line[4] == 'Y' || (line[4] >= '1' && line[4] <= '9')) {
                    string chrom = line.substr(1, l - 1) ;
                    //cout << "Collecting " << chrom << ".." << endl ;
                    while(std::getline(fasta_file, line)) {
                         if (line[0] == '>') {
                             break ;
                         }
                         for (int i = 0; i < line.length(); i++) {
                            line[i] = toupper(line[i]) ;
                         }
                         memcpy(buffer + n, line.c_str(), line.length()) ;
                         n += line.length() ;
                    }
                    buffer[n] = '\0' ;
                    cout << "Extracted " << std::setw(6) << std::left << chrom << " with " << std::setw(11) << left << n << " bases." << endl ;
                    char* s = (char*) malloc(sizeof(char) * (n + 1)) ;
                    memcpy(s, buffer, n + 1) ;
                    chromosomes.push_back(chrom) ;
                    chromosome_seqs[chrom] = s ;
                    n = 0 ;
                    continue ;
                }
        } 
        if (!std::getline(fasta_file, line)) {
            break ;
        }
    }
    free(buffer) ;
}
