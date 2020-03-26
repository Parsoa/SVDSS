#include <ctime>
#include <chrono>
#include <locale>
#include <string>
#include <thread>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <algorithm>
#include <unordered_map>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "sdsl/suffix_trees.hpp"

#define DEBUG 0

using namespace std ;
using namespace sdsl ;

typedef cst_sct3<csa_bitcompressed<int_alphabet<> > > cst_t ;
typedef cst_t::string_type string_type ;
typedef cst_t::char_type char_type ;
typedef std::vector<uint16_t> read_type ;
typedef uint16_t base_type ;

// prototypes

void traverse_cst(cst_t* cst) ;
void calculate_child_diff(cst_t* father, cst_t* mother, std::string child, std::vector<read_type*>* source) ;
const char* intergralize_string(read_type* source) ;
int search_sequence(cst_t* cst, read_type* seq, read_type::iterator begin, read_type::iterator end) ;
int search_sequence_backward(cst_t* cst, read_type* seq, read_type::iterator begin, read_type::iterator end) ;
//bool search_sequence(cst_t* cst, read_type* seq) ;

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
            continue ;
        }
        if (alignment->core.tid == -1) {
            continue ;
        }
        char* contig = bam_header->target_name[alignment->core.tid] ; //contig name (chromosome)
        if (strcmp(contig, "chr22") != 0) {
            continue ;
        }
        reads->push_back(new read_type()) ;
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < l; i++) {
            uint16_t base = uint16_t(seq_nt16_str[bam_seqi(q, i)]) ;
            (*reads)[u]->push_back(base) ;
            if (i == 20) {
                break ;
            }
        }
        (*reads)[u]->push_back(u + uint16_t(255 + 1)) ;
        n += 1 ;
        u += 1 ;
        if (n == 10) {
            break ;
            n = 0 ;
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

std::vector<read_type*>* process_bam_debug_child(string bam) {
    std::vector<read_type*>* reads = new std::vector<read_type*>() ;
    reads->push_back(new read_type({71, 67, 67, 84, 84, 67, 84, 67, 84, 84, 67, 65, 84, 71, 71, 65, 71, 67, 84, 67, 67, 256})) ;
    return reads ;
}

std::vector<read_type*>* process_bam_debug_father(string bam) {
    std::vector<read_type*>* reads = new std::vector<read_type*>() ;
    reads->push_back(new read_type({84, 67, 65, 84, 71, 71, 65, 71, 67, 84, 67, 67, 84, 84, 67, 84, 84, 71, 71, 67, 84, 67, 67, 71, 65, 71, 67, 71, 67, 67, 256})) ;
    return reads ;
}

cst_t* create_suffix_tree(std::string sample) {
    std::vector<read_type*>* reads = process_bam_debug_father(sample) ;
    int i = 0 ;
    int l = 0 ;
    int n = 0 ;
    int s = 0 ;
    std::vector<char*>* int_reads = new std::vector<char*>() ;
    for (auto it = reads->begin(); it != reads->end(); it++) { //iterate over reads
        std::vector<string> tmp ;
        l = 0 ;
        for (auto itt = (*it)->begin(); itt != (*it)->end(); itt++) { //iterate over characters in each read
            std::string b = std::to_string(*itt) ;
            tmp.push_back(b) ;
            l += b.length() + 1 ; // add one for space
        }
        i = 0 ;
        char* read = (char*) malloc(sizeof(char) * (l + 1)) ;
        for (auto itt = tmp.begin(); itt != tmp.end(); itt++) {
            strncpy(read + i, (*itt).c_str(), (*itt).length()) ; 
            read[i + (*itt).length()] = ' ' ;
            i += (*itt).length() + 1 ;
        }
        read[l] = '\0' ;
        cout << "|" << read << "|" << endl ;
        n += 1 ;
        int_reads->push_back(read) ;
    }
    cout << "Assembling master read.." << endl ;
    s = 0 ;
    for (auto it = int_reads->begin(); it != int_reads->end(); it++) {
        s += strlen(*it) ;
    }
    cout << "Master read size: " << s << " bytes.." << endl ;
    char* master_read = (char*) malloc(sizeof(char) * (s)) ; // NUll terminator overwrites last space
    i = 0 ;
    for (auto it = int_reads->begin(); it != int_reads->end(); it++) {
        strcpy(master_read + i, *it) ;
        i += strlen(*it) ;
    }
    master_read[s - 1] = '\0' ;
    cout << "Assembling tree.." << endl ;
    cst_t* cst = new cst_t() ;
    construct_im(*cst, (const char*)master_read, 'd') ;
    //csXprintf(cout, "%2I %3S %:4T", *cst);
    return cst ;
}

// from https://github.com/simongog/sdsl-lite/blob/master/tutorial/cst-traversal.cpp
void traverse_cst(cst_t* cst) {
    for (auto it = cst->begin(); it != cst->end(); ++it) {
        if (it.visit() == 1) { // node visited the first time
            auto v = *it ; // get the node by dereferencing the iterator
            if (cst->depth(v) <= 1000) { // if depth node is <= max_depth
                // process node, e.g. output it in format d-[lb, rb]
                cout << cst->depth(v) << "-[" << cst->lb(v) << "," << cst->rb(v) << "]" << endl ;
            } else { // skip the subtree otherwise
                it.skip_subtree() ;
            }
        }
    }
}

const char* intergralize_string(read_type* source) {
    int l = 0 ;
    int i = 0 ;
    std::vector<string> tmp ;
    for (auto it = (source)->begin(); it != (source)->end(); it++) { //iterate over characters in each read
        std::string b = std::to_string(*it) ;
        tmp.push_back(b) ;
        l += b.length() + 1 ; // add one for space
    }
    char* read = (char*) malloc(sizeof(char) * l) ;
    for (auto itt = tmp.begin(); itt != tmp.end(); itt++) {
        strncpy(read + i, (*itt).c_str(), (*itt).length()) ; 
        read[i + (*itt).length()] = ' ' ;
        i += (*itt).length() + 1 ;
    }
    read[l - 1] = '\0' ; //overwrites the final space
    return (const char*) read ;
}

void calculate_child_diff(cst_t* father, cst_t* mother, std::string child) {
    std::vector<read_type*>* reads = process_bam_debug_child(child) ;
    std::vector<char*> diff ;
    int n = reads->size() ;
    int m = 0 ;
    std::vector<read_type*>* mismatches = new std::vector<read_type*>() ;
    for (auto it = reads->begin(); it != reads->end(); it++) {
        cout << "------- matching -------" << endl ;
        int q = 0 ;
        int offset = 0 ;
        int l = (*it)->size() - 1 ; // ignore the terminator
        cout << "read length " << l << endl ;
        while (true) {
            std::this_thread::sleep_for (std::chrono::seconds(1));
            offset = search_sequence_backward(father, *it, (*it)->begin(), (*it)->end() - q - 1 - 1) ; // end() is one past the terminator, subtract two to get to the last base pair
            if (offset != -1) {
                int p = 0 ;
                cout << "binary search for longest mismatch at offset " << offset << ", " << (*it)->at(l - q - offset) << endl ;
                while (true) {
                    std::this_thread::sleep_for (std::chrono::seconds(1));
                    p += 1 ;
                    int end_offset = - q + p ;
                    int begin_offset = - q - offset - p ;
                    int end = l - 1 + end_offset <= l - 1 ? l - 1 + end_offset : l - 1 ;
                    int begin = l - 1 + begin_offset >= 0 ? l - 1 + begin_offset : 0 ;
                    cout << "begin offset: " << begin_offset << " , end offset: " << end_offset << endl ;
                    cout << "interval [" << begin << ", " << end << "]" << endl ;
                    int m = search_sequence(father, *it, (*it)->begin() + begin, (*it)->begin() + end) ;
                    if (m == 0) {
                        // check if subsequences of another string
                        cout << "matched." << endl ;
                        p -= 1 ;
                        break ;
                    }
                    if (begin == 0) {
                        break ;
                    }
                }
                q += offset + p ;
                if (q == l - 1) { // consumed whole read
                    break ;
                }
            } else {
                // matches whole read
                break ;
            }
        }
    }
    cout << m << " out of " << n << " not found." << endl ;
}

int search_sequence_backward(cst_t* cst, read_type* seq, read_type::iterator begin, read_type::iterator end) {
    uint64_t lb = 0, rb = cst->size() - 1 ;
    const char* read = intergralize_string(seq) ;
    read_type* tmp = new read_type(begin, end + 1) ;
    cout << "searching for: " << intergralize_string(tmp) << endl ;
    int offset = 0 ;
    // TODO: doesn't match last base
    for (auto itt = end; itt != begin and lb <= rb;) {
        offset++ ;
        if (backward_search(cst->csa, lb, rb, (char_type)*itt, lb, rb) > 0) {
            cout << *itt << endl ;
            --itt ;
            continue ;
        } else {
            cout << "Not found. offset " << offset << ", value " << *itt << endl ;
            return offset ;
        } 
    }
    //cout << "Found." << endl ;
    return -1 ;
} 

int search_sequence(cst_t* cst, read_type* seq, read_type::iterator begin, read_type::iterator end) {
    uint64_t lb = 0, rb = cst->size() - 1 ;
    // slice read
    read_type* tmp = new read_type(begin, end + 1) ;
    // convert to char array
    const char* read = intergralize_string(tmp) ;
    cout << "searching for: " << read << endl ;
    // copy to a vector for iterator access
    std::vector<char> c ;
    c.assign(read, read + strlen(read)) ;
    backward_search(cst->csa, lb, rb, c.begin(), c.end(), lb, rb) ;
    int match_size = rb + 1 - lb - tmp->size() ;
    cout << rb << ", " << lb << ", " << tmp->size() << endl ;
    return match_size ;
}

int main(int argc, char** argv) {
    string child = argv[3] ;
    string father = argv[1] ; 
    string mother = argv[2] ;
    cout << "Creating father suffix tree.." << endl ;
    cst_t* father_cst = create_suffix_tree(father) ;
    cout << "Creating mother suffix tree.." << endl ;
    //cst_t* mother_cst = create_suffix_tree(mother) ;
    calculate_child_diff(father_cst, father_cst, child) ;
    calculate_child_diff(father_cst, father_cst, father) ;
}