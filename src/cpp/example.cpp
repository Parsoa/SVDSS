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

#include "sdsl/suffix_trees.hpp"

#include "json.hpp"

#ifdef DEBUG_BUILD
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

#define NUM_THREADS 1

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

cst_t* create_suffix_tree(std::string sample) {
    //return nullptr ;
    std::vector<read_type*>* reads = process_bam(sample) ;
    int i = 0 ;
    int l = 0 ;
    int n = 0 ;
    int s = 0 ;
    cout << "Assembling master read.." << endl ;
    read_type* master_read = new read_type() ;
    for (auto it = reads->begin(); it != reads->end(); it++) { //iterate over reads
        l = 0 ;
        for (auto itt = (*it)->begin(); itt != (*it)->end(); itt++) { //iterate over characters in each read
            master_read->push_back(*itt) ;
        }
    }
    cout << "Master read size: " << master_read->size() << " bytes.." << endl ;
    cout << "Assembling tree.." << endl ;
    cst_t* cst = new cst_t() ;
    construct_im(*cst, master_read, 'd') ;
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
    NEBUG(cout_mutex.lock() ;)
    nlohmann::json payload ;
    cout << "dumping novel sequences..." << endl ;
    cout << "found " << seqs->size() << " novel sequences in the child." << endl ;
    string p = "batch_" + std::to_string(index) + ".json" ;
    std::ofstream o(p);
    nlohmann::json j(*seqs) ;
    o << j.dump(4) << std::endl ;
    cout << "done" << endl ;
    NEBUG(cout_mutex.unlock() ;)
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
        if (n == 10) {
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
                cout << "\r" ;
                cout << std::left << "thread " << setw(3) << index ;
                cout << " processed " << setw(8) << u << " reads, " ;
                cout << " took: " << setw(7) << std::fixed << s - t ;
                cout << " reads per second: " << u / (s - t) ;
                cout << " progress: " << setw(10) << std::fixed << p ;
                cout << " ETA: " << setw(12) << e ;
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
    read_type* read = new read_type({71, 67, 67, 84, 84, 67, 84, 67, 84, 84, 67, 65, 84, 71, 71, 65, 71, 65, 84, 67, 67, 84, 67, 65, 71, 67, 84, 65, 84, 71, 65, 256}) ;
    cst_t* cst = new cst_t() ;
    construct_im(*cst, (const char*) master_read, 'd') ;
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
