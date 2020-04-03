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
#include <mutex>
#include <unordered_map>
#include <pthread.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "sdsl/suffix_trees.hpp"
#include "sdsl/io.hpp"

//#include "suffixtree.h"

#include "json.hpp"

#ifdef DEBUG_BUILD
#  define DEBUG(x) x
#  define NEBUG(x)
#  define NUM_THREADS 1
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#  define NUM_THREADS 16
#endif


using namespace std ;
using namespace sdsl ;

// typedefs

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

// gobals
std::mutex cout_mutex ;

std::vector<read_type*>* process_bam(string bam) {
    samFile *bam_file = hts_open(bam.c_str(), "r") ;
    bam_hdr_t *bam_header = sam_hdr_read(bam_file) ; //read header
    bam1_t *alignment = bam_init1(); //initialize an alignment
    int n = 0 ;
    uint16_t u = 0 ;
    uint32_t len = 0 ;
    char* line ;
    // Timing
    time_t t ;
    time(&t) ;
    cout << "Processig BAM file" << endl ;
    std::vector<read_type*>* reads = new std::vector<read_type*>() ;
    //return reads ;
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
            DEBUG(if (i == 30) { break ; })
            //if (i == 30) { break ;}
        }
        (*reads)[u]->push_back(u + uint16_t(255 + 1)) ;
        n += 1 ;
        u += 1 ;
        DEBUG(if (n == 10) { break ; })
        if (n == 100) {
            //break ;
            n = 0 ;
            time_t s ;
            time(&s) ;
            cout.precision(10) ;
            if (s - t != 0) {
                cout << std::left << "processed " << setw(12) << u << " reads, " ;
                cout << " took: " << setw(7) << std::fixed << s - t ;
                cout << " reads per second: " << u / (s - t) << "\r" ;
            }
        }
    }
    cout << endl ;
    bam_destroy1(alignment) ;
    sam_close(bam_file) ;
    cout << "Extracted " << reads->size() << " reads.." << std::endl ;
    return reads ;
}

//cst_t* create_suffix_tree(std::string sample) {
//    //return nullptr ;
//    std::vector<read_type*>* reads = process_bam(sample) ;
//    cout << "Assembling master read.." << endl ;
//    read_type master_read ; //= new read_type() ;
//    for (auto it = reads->begin(); it != reads->end(); it++) { //iterate over reads
//        for (auto itt = (*it)->begin(); itt != (*it)->end(); itt++) { //iterate over characters in each read
//            master_read.push_back(*itt) ;
//        }
//    }
//    cout << "Master read size: " << master_read.size() << " bytes.." << endl ;
//    cout << "Assembling tree.." << endl ;
//    cst_t* cst = new cst_t() ;
//    std::string tmp_file = "master_read.txt" ;
//    store_to_file(master_read, tmp_file);
//    construct_im(*cst, master_read, 0) ;
//    return cst ;
//}

cst_t* create_suffix_tree(std::string sample) {
    //return nullptr ;
    std::vector<read_type*>* reads = process_bam(sample) ;
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
        DEBUG(cout << "|" << read << "|" << endl ;)
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

void print_read(read_type* read, read_type::iterator begin, read_type::iterator end) {
    for (auto t = begin; t != end; t++) {
        cout << *t << " " ;
    }
    cout << endl ;
}

void print_read(read_type* read) {
    print_read(read, read->begin(), read->end()) ;
}

// TODO: optimize this
std::string hash_string(read_type* source, read_type::iterator begin, read_type::iterator end) {
    int l = 0 ;
    int i = 0 ;
    std::vector<string> tmp ;
    for (auto it = begin; it != end; it++) { //iterate over characters in each read
        std::string b = std::to_string(*it) ;
        tmp.push_back(b) ;
        l += b.length() + 1 ; // add one for space
    }
    char read[l] ;
    for (auto itt = tmp.begin(); itt != tmp.end(); itt++) {
        strncpy(read + i, (*itt).c_str(), (*itt).length()) ;
        read[i + (*itt).length()] = ' ' ;
        i += (*itt).length() + 1 ;
    }
    read[l - 1] = '\0' ; //overwrites the final space
    string s = std::string(read) ;
    DEBUG(cout << "hash: " << s << endl ;)
    return s ;
}

std::string hash_string(read_type* source) {
    return hash_string(source, source->begin(), source->end()) ;
}

void output_diff(std::string path, int index, std::unordered_map<std::string, int>* seqs) {
    //NEBUG(cout_mutex.lock() ;)
    nlohmann::json payload ;
    cout << "dumping novel sequences..." << endl ;
    cout << "found " << seqs->size() << " novel sequences in the child." << endl ;
    string p = "batch_" + std::to_string(index) + ".json" ;
    std::ofstream o(p);
    nlohmann::json j(*seqs) ;
    o << j.dump(4) << std::endl ;
    cout << "done" << endl ;
    //NEBUG(cout_mutex.unlock() ;)
}

struct thread_data {
    int index ;
    cst_t* father ;
    cst_t* mother ;
    std::string child ;
    std::vector<read_type*>* reads ;
} ;

void* calculate_child_diff_t(void* args) {
    struct thread_data* t_data = (struct thread_data*) args ;
    int index = t_data->index ;
    cst_t* father = t_data->father ;
    cst_t* mother = t_data->mother ;
    //cout << "starting diff thread " << index << ".." << endl ;
    time_t t ;
    time(&t) ;
    int n = 0 ;
    int u = 0 ;
    std::vector<read_type*>* reads = t_data->reads ;
    std::unordered_map<std::string, int>* mismatched_strings = new std::unordered_map<std::string, int>() ;
    for (auto it = reads->begin(); it != reads->end(); it++) {
    //while(true) {
        //auto it = reads->begin() ;
        u += 1 ;
        DEBUG(cout << "------- matching -------" << endl ;)
        if (u % NUM_THREADS != index) {
            continue ;
        }
        int q = 0 ;
        int offset = 0 ;
        int l = (*it)->size() - 1 ; // ignore the terminator
        DEBUG(cout << "read length " << l << endl ;)
        while (true) {
            DEBUG(std::this_thread::sleep_for (std::chrono::seconds(1));)
            DEBUG(cout << "continue at offset " << q << endl ;)
            offset = search_sequence_backward(father, *it, (*it)->begin(), (*it)->end() - 1 - 1 - q + 1) ; // end() is one past the terminator, subtract two to get to the last base pair
            if (offset != -1) {
                DEBUG(cout << "binary search for longest mismatch at offset " << offset << ", " << (*it)->at(l - q - offset) << endl ;)
                // (end of read) - (offset into the read) - (change since last offset) + (adjustment) 
                int pivot = (l - 1) - q - offset + 1 ;
                int end_limit = pivot;
                int begin_limit = pivot ;
                //
                int end = pivot ;
                int begin = pivot ;
                // fix the end
                while (end <= l - 1) {
                    begin = pivot ;
                    DEBUG(std::this_thread::sleep_for (std::chrono::seconds(1));)
                    DEBUG(cout << "interval [" << begin << ", " << end << "]" << endl ;)
                    int m = search_sequence(father, *it, (*it)->begin() + begin, (*it)->begin() + end + 1) ;
                    if (m == 0) {
                        DEBUG(cout << "added" << endl ;)
                        if (end > pivot) {
                            end_limit = end ;
                            mismatched_strings->emplace(std::make_pair(hash_string(*it, (*it)->begin() + begin, (*it)->begin() + end + 1), end - begin + 1)) ;
                            break ;
                        }
                    }
                    end += 1 ;
                }
                while (end > pivot) {
                    end -= 1 ;
                    begin = begin_limit ;
                    while (begin >= 0) {
                        DEBUG(std::this_thread::sleep_for (std::chrono::seconds(1));)
                        DEBUG(cout << "interval [" << begin << ", " << end << "]" << endl ;)
                        int m = search_sequence(father, *it, (*it)->begin() + begin, (*it)->begin() + end + 1) ;
                        if (m == 0) {
                            DEBUG(cout << "added" << endl ;)
                            if (begin < pivot) {
                                begin_limit = begin ;
                                mismatched_strings->emplace(std::make_pair(hash_string(*it, (*it)->begin() + begin, (*it)->begin() + end + 1), end - begin + 1)) ;
                                break ;
                            }
                        }
                        begin -= 1 ;
                    }
                }
            } else {
                break ;
            }
            q += offset ;
            if (q >= l - 1) {
                break ;
            }
        }
        n += 1 ;
        if (n == 100) {
            n = 0 ;
            time_t s ;
            time(&s) ;
            cout.precision(6) ;
            double p = u / (double(reads->size()) / NUM_THREADS) ;
            double e = (((1.0 - p) * (s - t)) / p) / 3600 ;
            if (s - t != 0) {
                //NEBUG(cout_mutex.lock() ;)
                for (int j = 0; j < NUM_THREADS - index; j++) {
                    cout << "\x1b[A" ;
                }
                cout << "\r"
                << std::left << "thread " << setw(3) << index
                << " processed " << setw(8) << u << " reads, "
                << " took: " << setw(7) << std::fixed << s - t
                << " reads per second: " << u / (s - t)
                << " progress: " << setw(10) << std::fixed << p
                << " ETA: " << setw(12) << e ;
                for (int j = 0; j < NUM_THREADS - index; j++) {
                    cout << endl ;
                }
                //NEBUG(cout_mutex.unlock() ;)
            }
        }
    }
    output_diff("", index, mismatched_strings) ;
    pthread_exit(NULL) ;
}

void calculate_child_diff(cst_t* father, cst_t* mother, std::string child) {
    pthread_t threads[NUM_THREADS] ;
    struct thread_data* t_data[NUM_THREADS] ;
    std::vector<read_type*>* reads = process_bam(child) ;
    cout << "======================== STATUS ============================ " << endl ;
    for (int t = 0; t < NUM_THREADS; t++) {
        cout << endl ;
    }
    cout << endl ;
    for (int t = 0; t < NUM_THREADS; t++) {
        t_data[t] = new thread_data() ;
        t_data[t]->index = t ;
        t_data[t]->father = father ;
        t_data[t]->mother = mother ;
        t_data[t]->child = child ;
        t_data[t]->reads = reads ;
        int rc = pthread_create(&threads[t], NULL, calculate_child_diff_t, (void *) t_data[t]) ;
    }
    int rc ;
    void *status ;
    for (int t = 0; t < NUM_THREADS; t++) {
        rc = pthread_join(threads[t], &status);
    }
}

// The end iterator passed to this should +1 the last desired base pair
int search_sequence_backward(cst_t* cst, read_type* seq, read_type::iterator begin, read_type::iterator end) {
    uint64_t lb = 0, rb = cst->size() - 1 ;
    DEBUG(cout << "searching for: " ;
    for (auto t = begin; t != end; t++) {
        cout << *t << " " ;
    }
    cout << endl ;)
    int offset = 0 ;
    if (begin != seq->begin()) {
        begin-- ;
    }
    end-- ;
    for (auto itt = end; itt != begin and lb <= rb;) {
        offset++ ;
        if (backward_search(cst->csa, lb, rb, (char_type)*itt, lb, rb) > 0) {
            DEBUG(cout << *itt << endl ;)
            --itt ;
            continue ;
        } else {
            DEBUG(cout << "Not found. offset " << offset << ", value " << *itt << endl ;)
            return offset ;
        } 
    }
    return -1 ;
} 

// The end iterator passed to this should +1 the last desired base pair
int search_sequence(cst_t* cst, read_type* seq, read_type::iterator begin, read_type::iterator end) {
    uint64_t lb = 0, rb = cst->size() - 1 ;
    DEBUG(cout << "searching for: " ;
    for (auto t = begin; t != end; t++) {
        cout << *t << " " ;
    }
    cout << endl ;)
    backward_search(cst->csa, lb, rb, begin, end, lb, rb) ;
    int match_num = rb + 1 - lb ;
    return match_num ;
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
    //calculate_child_diff(father_cst, father_cst, father) ;
}
