#include <locale>
#include <ctime>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <cstdint>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <bitset>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "sdsl/suffix_trees.hpp"

#define DEBUG 0

using namespace std ;
using namespace sdsl ;

//typedef cst_sct3<> cst_t ;
typedef cst_sct3<csa_bitcompressed<int_alphabet<> > > cst_t ;
typedef cst_t::string_type string_type ;
typedef cst_t::char_type char_type ;
typedef std::vector<uint16_t> read_type ;
typedef uint16_t base_type ;

// prototypes

void calculate_child_diff(cst_t* father, cst_t* mother, std::string child, std::vector<read_type>* source) ;

//

std::vector<read_type*>* process_bam(string bam) {
    samFile *bam_file = hts_open(bam.c_str(), "r") ;
    bam_hdr_t *bam_header = sam_hdr_read(bam_file) ; //read header
    bam1_t *alignment = bam_init1(); //initialize an alignment
    int n = 0 ;
    uint16_t u = 0 ;
    uint32_t len = 0 ;
    char* line ; //= (char*) malloc(200) ;
    // Timing
    time_t t ;
    time(&t) ;
    cout << "Processig BAM file" << endl ;
    std::vector<read_type*>* reads = new std::vector<read_type*>() ;
    while (sam_read1(bam_file, bam_header, alignment) > 0){
        uint32_t l = alignment->core.l_qseq ; //length of the read
        char qname[alignment->core.l_qname] ;
        strncpy(qname, (char*) alignment->data, alignment->core.l_qname) ;
        if (l <= 1) {
            cout << qname << " " << alignment->core.tid << endl ;
            continue ;
        }
        if (alignment->core.tid == -1) {
            continue ;
        }
        char* contig = bam_header->target_name[alignment->core.tid] ; //contig name (chromosome)
        //cout << " " <<  contig << " " << "\r" ;
        if (strcmp(contig, "chr22") != 0) {
            continue ;
        }
        //line = (uint16_t*) malloc(l + 2) ;
        reads->push_back(new read_type()) ;
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < l; i++) {
            (*reads)[u]->push_back(uint16_t(seq_nt16_str[bam_seqi(q, i)])) ;
            //line[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        //reads[u].push_back(u + 255 + 1) ;
        (*((*reads)[u]))[20] = u + uint16_t(255 + 1) ;
        //line[l] = '$' ;
        //line[l + 1] = '\0' ;
        //reads->push_back(line) ;
        n += 1 ;
        u += 1 ;
        if (n == 1000) {
            n = 0 ;
            break ;
            time_t s ;
            time(&s) ;
            cout.precision(10) ;
            if (s - t != 0 && DEBUG == 0) {
                cout << std::left << "processed " << setw(12) << u << " reads, " ;
                cout << " took: " << setw(7) << std::fixed << s - t ;
                cout << " reads per second: " << u / (s - t) << "\r" ;
            }
        }
    }
    cout << endl ;
    bam_destroy1(alignment) ;
    sam_close(bam_file) ;
    cout << "Extracted " << u << " reads.." << std::endl ;
    return reads ;
}

cst_t* create_suffix_tree(std::string sample) {
    //std::vector<char*>* reads = process_bam(sample) ;
    std::vector<read_type*>* reads = process_bam(sample) ;
    int s = 0 ;
    for (auto it = reads->begin(); it != reads->end(); it++) {
        //s += strlen(*it) ;
        s += 20 ;
    }
    int i = 0 ;
    int n = 0 ;
    char* a = (char*) malloc(sizeof(char) * (s + 1)) ;
    for (auto it = reads->begin(); it != reads->end(); it++) {
        n += 1 ;
        cout << "Adding " << n << "th read.." << "\r" ;
        //int l = strlen(*it) ;
        int l = 20 ;
        memcpy(a + i, *it, l * sizeof(uint16_t)) ;
        i += l ;
    }
    cout << endl ;
    cout << "Assembling tree.." << endl ;
    a[s] = '$' ;
    a[s] = '\0' ;
    cst_t* cst = new cst_t() ;
    construct_im(*cst, (const char*)a, 1) ;
    //calculate_child_diff(cst, cst, sample, reads) ;
    return cst ;
}

void calculate_child_diff(cst_t* father, cst_t* mother, std::string child, std::vector<read_type*>* source) {
    std::vector<read_type*>* reads = process_bam(child) ;
    std::vector<char*> diff ;
    for (auto it = reads->begin(); it != reads->end(); it++) {
        uint64_t lb = 0, rb = father->size() - 1 ;
        /*char subbuff[20] ;
        memcpy(subbuff, *it, 19) ;
        subbuff[19] = '\0' ;
        cst_t::string_type pat = subbuff ;
        cout << strlen(*it) << endl ;
        cout << pat << endl ;
        backward_search(father->csa, lb, rb, pat.begin(), pat.end(), lb, rb);
        cout << "size = " << rb + 1 - lb << endl ;
        int b = 0 ;
        for (auto itt = source->begin(); itt != source->end(); itt++) {
            if (strcmp(*itt, *it) == 0) {
                b = 1 ;
                cout << "Found.." << endl ;
                break ;
            }
        }
        if (b == 0) {
            cout << "Not found.." << endl ;
        }
        cout << "-----" << endl ;
        cst_t::string_type pat = *it ;
        for (auto itt = pat.end(); itt != pat.begin() and lb <= rb;) {
            --itt ;
            if (backward_search(father->csa, lb, rb, (typename cst_t::char_type)*itt, lb, rb) > 0) {
                cout << "[" << lb << "," << rb << "]" << endl;
                cout << "matched " << *itt << endl;
            } else {
                cout << "Not matched " << *itt << endl;
            }
        }*/
    }
}

int main(int argc, char** argv) {
    /*string path(argv[2]) ;
    string fastq(argv[3]) ;
    int index = std::stoi(string(argv[1]), nullptr, 10) ;
    int threads = std::stoi(string(argv[4]), nullptr, 10) ;
    JOB = std::stoi(string(argv[5]), nullptr, 10) ;
    DEBUG = std::stoi(string(argv[6]), nullptr, 10) ; */
    string father = "/share/hormozdiarilab/Codes/Stella/data/PUR-HiFi/HG00731/HG00731_hifi_r54329U_20190531_173856_A01.bam.2" ;
    string mother = "/share/hormozdiarilab/Codes/Stella/data/PUR-HiFi/HG00732/HG00732_hifi_r54329U_20190604_223930_A01.sam.1" ;
    string child = "/share/hormozdiarilab/Codes/Stella/data/PUR-HiFi/HG00733/HG00733_hifi_r54329U_20190827_172128_1_A01.sam.1" ;
    cout << "Creating father suffix tree.." << endl ;
    cst_t* father_cst = create_suffix_tree(father) ;
    //cout << "Creating mother suffix tree.." << endl ;
    //cst_sct3<>* mother_cst = create_suffix_tree(mother) ;
    //calculate_child_diff(father_cst, father_cst, father) ;
}
