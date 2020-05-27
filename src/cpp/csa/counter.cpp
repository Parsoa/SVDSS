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
typedef std::vector<uint32_t> read_type ;
typedef uint32_t base_type ;

base_type terminus = 85 ;

std::string chrom ;
std::string output_path ;

// prototypes

void traverse_cst(cst_t* cst) ;
void calculate_child_diff(cst_t* father, cst_t* mother, std::string child, std::vector<read_type*>* source) ;
const char* intergralize_string(read_type* source) ;
int search_sequence(cst_t* cst, read_type* seq, read_type::iterator begin, read_type::iterator end) ;
int search_sequence_backward(cst_t* cst, read_type* seq, read_type::iterator begin, read_type::iterator end) ;

// gobals
std::mutex cout_mutex ;

//

void memory_usage() {
    double vm_usage = 0.0;
    double resident_set = 0.0;
    unsigned long vsize ;
    long rss ;
    std::string ignore ;
    std::ifstream ifs("/proc/self/stat", std::ios_base::in) ;
    ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> vsize >> rss ;
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024 ; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0 ;
    resident_set = rss * page_size_kb ;
    cout << "VM: " << vm_usage << "; RSS: " << resident_set << endl ;
}


std::vector<read_type*>* process_bam(string bam) {
    samFile *bam_file = hts_open(bam.c_str(), "r") ;
    bam_hdr_t *bam_header = sam_hdr_read(bam_file) ; //read header
    char **target_name = bam_header->target_name ;
    int32_t n_targets = bam_header->n_targets ;
    uint32_t *target_len = bam_header->target_len ;
    int32_t chr21 = 0 ;
    for (int i = 0 ; i < n_targets ; i++) {
        //cout << target_name[i] << endl ;
        if (strcmp(target_name[i], chrom.c_str()) == 0) {
            chr21 = i ;
            cout << "chr21: " << i << endl ;
        }
    }
    bam1_t *alignment = bam_init1(); //initialize an alignment
    int e = 0 ;
    int n = 0 ;
    int r = 0 ;
    uint32_t u = 0 ;
    uint32_t len = 0 ;
    uint64_t sum = 0 ;
    uint64_t sum_v = 0 ;
    char* line ;
    std::vector<read_type*>* reads = new std::vector<read_type*>() ;
    // Timing
    time_t t ;
    time(&t) ;
    cout << "Processig BAM file" << endl ;
    //return reads ;
    while (sam_read1(bam_file, bam_header, alignment) > 0){
        uint32_t l = alignment->core.l_qseq ; //length of the read
        char qname[alignment->core.l_qname] ;
        strncpy(qname, (char*) alignment->data, alignment->core.l_qname) ;
        r += 1 ;
        if (l <= 1) {
            e += 1 ;
            continue ;
        }
        if (alignment->core.tid == chr21) {
            sum += l ;
        } else {
            continue ;
        }
        //char* contig = bam_header->target_name[alignment->core.tid] ; //contig name (chromosome)
        //if (strcmp(contig, "chr21") != 0) {
        //    continue ;
        //}
        //sum += l ;
        reads->push_back(new read_type()) ;
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < l; i++) {
            uint32_t base = uint32_t(seq_nt16_str[bam_seqi(q, i)]) ;
            (*reads)[u]->push_back(base) ;
            DEBUG(if (i == 30) { break ; })
        }
        (*reads)[u]->push_back(terminus) ;
        terminus += 1 ;
        sum_v += (*reads)[u]->size() ;
        n += 1 ;
        u += 1 ;
        DEBUG(if (n == 10) { break ; })
        if (n == 100) {
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
    cout << "Processed " << r << " reads. Ignored " << e << " reads. Extracted " << reads->size() << " = " << u << " reads, totalling " << sum << " = " << sum_v << " bases. Terminus: " << terminus << std::endl ;
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
//

uint64_t encode_reads(std::vector<char*>* int_reads, std::vector<read_type*>* reads) {
    uint64_t l = 0 ;
    uint64_t s = 0 ;
    uint64_t r = 0 ;
    for (auto it = reads->begin(); it != reads->end(); it++) { //iterate over reads
        std::vector<string> tmp ;
        // calculate string size for read
        l = 0 ;
        r += (*it)->size() ;
        for (auto itt = (*it)->begin(); itt != (*it)->end(); itt++) { //iterate over characters in each read
            std::string b = std::to_string(*itt) ;
            tmp.push_back(b) ;
            l += b.length() + 1 ; // add one for space
        }
        uint64_t i = 0 ;
        char* read = (char*) malloc(sizeof(char) * l) ;
        //cout << "int size: " << (*it)->size() << ", char size: " << l + 1 << endl ;
        for (auto itt = tmp.begin(); itt != tmp.end(); itt++) {
            strncpy(read + i, (*itt).c_str(), (*itt).length()) ; 
            read[i + (*itt).length()] = ' ' ;
            i += (*itt).length() + 1 ;
        }
        //cout << i << " = " << l << endl ;
        read[l - 1] = '\0' ;
        DEBUG(cout << "|" << read << "|" << endl ;)
        s += l ;
        int_reads->push_back(read) ;
    }
    return s ;
}

//std::unordered_map<std::string, int> decode_reads(std::unordered_map<std::string, int>* reads) {
//    std::unordered_map<std::string, int> decoded_reads ;
//    for (auto it = reads->begin(); it != reads->end(); it++) { //iterate over reads
//        std::vector<string> tmp ;
//        // calculate string size for read
//        char* token = strtok(it.first, " ") ;
//        while (token != nullptr) {
//            if (strcmp(token, "65") == 0) {
//                tmp += "A" ;
//            } else if(strcmp(token, "67") == 0) {
//                tmp += "C" ;
//            } else if(strcmp(token, "71") == 0) {
//                tmp += "G" ;
//            } else if(strcmp(token, "84") == 0) {
//                tmp += "T" ;
//            } else {
//                break ;
//            }
//            char* token = strtok(NULL, " ") ;
//        }
//        decoded_reads[tmp] = it.second ;
//    }
//    return decode_reads ;
//}

//cst_t* create_suffix_tree(std::string father, std::string mother) {
//    std::vector<read_type*>* reads = process_bam(father) ;
//    uint64_t l = 0 ;
//    uint64_t s = 0 ;
//    uint64_t r = 0 ;
//    std::vector<char*>* int_reads = new std::vector<char*>() ;
//    for (auto it = reads->begin(); it != reads->end(); it++) { //iterate over reads
//        std::vector<string> tmp ;
//        // calculate string size for read
//        l = 0 ;
//        r += (*it)->size() ;
//        for (auto itt = (*it)->begin(); itt != (*it)->end(); itt++) { //iterate over characters in each read
//            std::string b = std::to_string(*itt) ;
//            tmp.push_back(b) ;
//            l += b.length() + 1 ; // add one for space
//        }
//        uint64_t i = 0 ;
//        char* read = (char*) malloc(sizeof(char) * (l + 1)) ;
//        //cout << "int size: " << (*it)->size() << ", char size: " << l + 1 << endl ;
//        for (auto itt = tmp.begin(); itt != tmp.end(); itt++) {
//            strncpy(read + i, (*itt).c_str(), (*itt).length()) ; 
//            read[i + (*itt).length()] = ' ' ;
//            i += (*itt).length() + 1 ;
//        }
//        //cout << i << " = " << l << endl ;
//        read[i] = '\0' ;
//        DEBUG(cout << "|" << read << "|" << endl ;)
//        s += l ;
//        int_reads->push_back(read) ;
//    }
//    cout << "Total number of bases: " << r << endl ;
//    cout << "Estimated master read size: " << s << " bytes.." << endl ;
//    cout << "Assembling master read.." << endl ;
//    char* master_read = (char*) malloc(sizeof(char) * (s)) ; // NUll terminator overwrites last space
//    uint64_t i = 0 ;
//    for (auto it = int_reads->begin(); it != int_reads->end(); it++) {
//        strcpy(master_read + i, *it) ;
//        i += strlen(*it) ;
//    }
//    std::ofstream of("master_read.txt") ;
//    cout << "Dumping master sequences..." << endl ;
//    of << master_read ;
//    of.close() ;
//    master_read[s - 1] = '\0' ; // NULL terminator overwrites the last space
//    cout << "Assembling tree.." << endl ;
//    cst_t* cst = new cst_t() ;
//    construct_im(*cst, (const char*)master_read, 'd') ;
//    return cst ;
//}

cst_t* create_suffix_tree(std::string father, std::string mother) {
    std::vector<read_type*>* reads_father = process_bam(father) ;
    std::vector<read_type*>* reads_mother = process_bam(mother) ;
    std::vector<char*>* int_reads = new std::vector<char*>() ;
    uint64_t s = 0 ;
    s += encode_reads(int_reads, reads_father) ;
    s += encode_reads(int_reads, reads_mother) ;
    cout << "Estimated master read size: " << s << " bytes.." << endl ;
    cout << "Assembling master read.." << endl ;
    char* master_read = (char*) malloc(sizeof(char) * (s)) ; // NUll terminator overwrites last space
    uint64_t i = 0 ;
    for (auto it = int_reads->begin(); it != int_reads->end(); it++) {
        strcpy(master_read + i, *it) ;
        i += strlen(*it) ;
    }
    master_read[s - 1] = '\0' ; // NULL terminator overwrites the last space
    memory_usage() ;
    cout << "Assembling tree.." << endl ;
    cst_t* cst = new cst_t() ;
    construct_im(*cst, (const char*)master_read, 'd') ;
    memory_usage() ;
    cout << "Freeing memory.." << endl ;
    delete int_reads ;
    delete reads_father ;
    delete reads_mother ;
    free(master_read) ;
    memory_usage() ;
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
// The end iterator passed to this should be +1 the last desired base pair
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
    string p = output_path + "/batch_" + std::to_string(index) + ".json" ;
    std::ofstream o(p);
    nlohmann::json j(*seqs) ;
    o << j.dump(4) << std::endl ;
    cout << "done" << endl ;
    //NEBUG(cout_mutex.unlock() ;)
}

struct thread_data {
    int index ;
    cst_t* cst ;
    std::string child ;
    std::vector<read_type*>* reads ;
} ;

void* calculate_child_diff_t(void* args) {
    struct thread_data* t_data = (struct thread_data*) args ;
    int index = t_data->index ;
    cst_t* cst = t_data->cst ;
    //cout << "starting diff thread " << index << ".." << endl ;
    time_t t ;
    time(&t) ;
    int n = 0 ;
    int u = 0 ;
    std::vector<read_type*>* reads = t_data->reads ;
    std::unordered_map<std::string, int>* mismatched_strings = new std::unordered_map<std::string, int>() ;
    for (auto it = reads->begin(); it != reads->end(); it++) {
        u += 1 ;
        DEBUG(cout << "------- matching -------" << endl ;)
        if (u % NUM_THREADS != index) {
            continue ;
        }
        int offset = 0 ;
        int l = (*it)->size() - 1 ; // ignore the terminator
        DEBUG(cout << "read length " << l << endl ;)
        offset = search_sequence_backward(cst, *it, (*it)->begin(), (*it)->end() - 1) ; // end() is one past the terminator, subtract two to get to the last base pair
        if (offset != -1) {
            DEBUG(std::this_thread::sleep_for (std::chrono::seconds(1));)
            //
            int end = l - 1 ;
            int begin = l - 1 - offset + 1 ;
            int prev_end = end ;
            while (begin >= 0) {
                DEBUG(cout << "== end: " << end << ", begin: " << begin << endl ;)
                // fix the end
                while (end > begin) {
                    DEBUG(std::this_thread::sleep_for (std::chrono::seconds(1));)
                    DEBUG(cout << "-E interval [" << begin << ", " << end << "]" << endl ;)
                    int m = search_sequence(cst, *it, (*it)->begin() + begin, (*it)->begin() + end + 1) ;
                    if (m == 0) {
                        prev_end = end ;
                    } else {
                        DEBUG(cout << "## fixed end at " << prev_end << endl ;)
                        mismatched_strings->emplace(std::make_pair(hash_string(*it, (*it)->begin() + begin, (*it)->begin() + prev_end + 1), prev_end - begin + 1)) ;
                        break ;
                    }
                    end -= 1 ;
                }
                begin -= 1 ;
                // fix the beginning
                while (begin >= 0) {
                    DEBUG(std::this_thread::sleep_for (std::chrono::seconds(1));)
                    DEBUG(cout << "-B interval [" << begin << ", " << end << "]" << endl ;)
                    int m = search_sequence(cst, *it, (*it)->begin() + begin, (*it)->begin() + end + 1) ;
                    if (m == 0) {
                        DEBUG(cout << "@@ fixed begin at " << begin << endl ;)
                        break ;
                    }
                    begin -= 1 ;
                }
                // repeat until we reach the beginning of the sequence
            }
        } else {
            break ;
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

void calculate_child_diff(cst_t* cst, std::string child) {
    pthread_t threads[NUM_THREADS] ;
    struct thread_data* t_data[NUM_THREADS] ;
    cout << "Loading child data..." << endl ;
    std::vector<read_type*>* reads = process_bam(child) ;
    cout << "======================== STATUS ============================ " << endl ;
    for (int t = 0; t < NUM_THREADS; t++) {
        cout << endl ;
    }
    cout << endl ;
    for (int t = 0; t < NUM_THREADS; t++) {
        t_data[t] = new thread_data() ;
        t_data[t]->cst = cst ;
        t_data[t]->index = t ;
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

// The end iterator passed to this should be +1 the last desired base pair
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

// The end iterator passed to this should be +1 the last desired base pair
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
    output_path = argv[4] ;
    chrom = argv[5] ;
    cout << "Creating father suffix tree.." << endl ;
    //cst_t* father_cst = create_suffix_tree(father) ;
    //cout << "Creating mother suffix tree.." << endl ;
    //cst_t* mother_cst = create_suffix_tree(mother) ;
    cst_t* parent_cst = create_suffix_tree(father, mother) ;
    calculate_child_diff(parent_cst, child) ;
    //calculate_child_diff(father_cst, father_cst, father) ;
}
