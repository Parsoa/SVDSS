#include <iostream>
#include <string>
#include <vector>

#include <zlib.h>

#include <omp.h>

#include "kseq.h"

#include "fermi.h"
#include "rld.h"

#include "kmc_api/kmc_file.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

extern "C" int main_build(int argc, char *argv[]);
extern "C" void seq_char2nt6(int l, unsigned char *s);

#ifdef NDEBUG
#define DEBUG(msg)
#else
#define DEBUG(msg) do { cerr << msg << endl; } while (0)
#endif

// FIXME: change string to char. Strings are more convenient for debugging
static const vector<string> int2char ({"$", "A", "C", "G", "T", "N"});

// TODO: use kseq instead of this (?)
struct fastq_entry {
  string head;
  string seq;
  string qual;
  uint start;
  uint len;

  fastq_entry(const string &h, const string &s, const string &q, const uint st = 0, const uint l = 0) : head(h), seq(s), qual(q) {
    start = st;
    len = l;
  }
};

// UTILITIES
// FIXME: this crashes if input is FASTA
fastq_entry get_solution(fastq_entry fqe, int s, int l) {
  string S (fqe.seq, s, l);
  string Q (fqe.qual, s, l);
  return fastq_entry(fqe.head, S, Q, s, l);
}

string interval2str(fmintv_t sai) {
  return "[" + to_string(sai.x[0]) + "," + to_string(sai.x[1]) + "," + to_string(sai.x[2]) + "]";
}

/**
 * Backward search a string in the index
 **/
bool backward_search(rld_t *index, const uint8_t *P, int p2) {
  fmintv_t sai; // fmintv_t is the struct used to store a SA interval.
  // typedef struct {
  //    uint64_t x[3]; // 0: start of the interval, backward; 1: forward; 2: size of the interval
  //    uint64_t info;
  // } fmintv_t;
  fm6_set_intv(index, P[p2], sai);
  while(sai.x[2] != 0 && p2 > 0) {
    --p2;
    fmintv_t osai[6];
    fm6_extend(index, &sai, osai, 1); //1: backward, 0: forward
    sai = osai[P[p2]];
  }
  return sai.x[2] != 0;
}

vector<fastq_entry> algo_f3(rld_t *index, fastq_entry fqe) {
  vector<fastq_entry> solutions;

  char *seq = new char[fqe.seq.size()+1];
  strcpy(seq, fqe.seq.c_str());
  uint8_t *P = (uint8_t*)seq;
  int m = fqe.seq.size();
  int p2 = m - 1;
  seq_char2nt6(m, P); // convert to integers

  fmintv_t sai;
  int bmatches = 0;
  fm6_set_intv(index, P[p2], sai);
  DEBUG("BS from " + int2char[P[p2]] + " (" + to_string(p2) + "): " + interval2str(sai));
  while(p2 >= 0) {
    while(sai.x[2] != 0 && p2 > 0) {
      ++bmatches;
      --p2;
      fmintv_t osai[6]; // output SA intervals (one for each symbol between 0 and 5)
      fm6_extend(index, &sai, osai, 1);
      sai = osai[P[p2]];
      DEBUG("'- BE with " + int2char[P[p2]] + " (" + to_string(p2) + "): " + interval2str(sai));
    }

    if(sai.x[2] != 0 && p2 <= 0) {
      DEBUG("E1.");
      break;
    }
    DEBUG("Mismatch " + int2char[P[p2]] + " (" + to_string(p2) + "). bmatches: " + to_string(bmatches));

    // FORWARD SEARCH
    int fp2 = p2;
    int fmatches = 0;
    fm6_set_intv(index, P[fp2], sai);
    fmintv_t prev_sai;
    DEBUG("FS from " + int2char[P[fp2]] + " (" + to_string(fp2) + "): " + interval2str(sai));
    while(sai.x[2] != 0 && fmatches < bmatches && fp2 < m-1) { // fmatches < bmatches && fp2 < m-1 can be useless
      ++fmatches;
      ++fp2;
      fmintv_t osai[6];
      fm6_extend(index, &sai, osai, 0);
      prev_sai = sai;
      sai = osai[P[fp2] >= 1 && P[fp2] <= 4 ? 5 - P[fp2] : P[fp2]];
      DEBUG("'- FE with " + int2char[P[fp2]] + " (" + to_string(fp2) + "): " + interval2str(sai));
    }
    DEBUG("Mismatch " + int2char[P[fp2]] + " (" + to_string(fp2) + "). fmatches: " + to_string(fmatches));

    int start_pos = fp2-fmatches; // here I can use p2
    int len = fmatches+1; // here I can use fp2-p2[+1]+1
    fastq_entry sol = get_solution(fqe, start_pos, len);
    solutions.push_back(sol);

    // here prev_sai should be initialized - in a real case scenario, I expect to forward search at least one char
    assert(prev_sai.x[2] != 0);
    // Note: I cannot fsearch till the end of the read
    sai = prev_sai;
    bmatches = fmatches;
    // --p2; // I will decrement it at the start of the next iteration
    if(p2<0)
      DEBUG("E2.");
    else
      DEBUG("CT " + int2char[P[p2]] + " (" + to_string(p2) + "): " + interval2str(sai));
  }

  delete[] seq;
  return solutions;
}

vector<vector<vector<fastq_entry>>> psearch_f3(const vector<fastq_entry> &entries, rld_t *index, int threads) {
  vector<vector<vector<fastq_entry>>> solutions (threads);

#pragma omp parallel for num_threads(threads)
  for(uint i=0; i<entries.size(); ++i) {
    vector<fastq_entry> sol = algo_f3(index, entries[i]);
    solutions[omp_get_thread_num()].push_back(sol);
  }

  return solutions;
}

/**
 * Find all strings satisfying formulation 3
 **/
int search_f3(int argc, char *argv[]) {
  char *index_path = argv[1];
  char *sample_path = argv[2];
  int threads = atoi(argv[3]);

  rld_t *index = rld_restore(index_path);

  gzFile fp = gzopen(sample_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  vector<fastq_entry> entries;
  while ((l = kseq_read(seq)) >= 0) {

    entries.push_back(fastq_entry(seq->name.s, seq->seq.s, seq->qual.s));
    if(entries.size() == 50000) {
      vector<vector<vector<fastq_entry>>> solutions = psearch_f3(entries, index, threads);

      for(const auto &th : solutions) {
	for(const auto sol : th) {
	  for(const auto s : sol) {
	    cout << "@" << s.head << ".css" << "_" << s.start << ":" << s.start+s.len-1 << endl
		 << s.seq << endl
		 << "+" << endl
		 << s.qual << endl;
	  }
	}
      }

      entries.clear();
    }
  }

  if(!entries.empty()) {
    vector<vector<vector<fastq_entry>>> solutions = psearch_f3(entries, index, threads);

    for(const auto &th : solutions) {
      for(const auto sol : th) {
	for(const auto s : sol) {
	  cout << "@" << s.head << ".css" << "_" << s.start << ":" << s.start+s.len-1 << endl
	       << s.seq << endl
	       << "+" << endl
	       << s.qual << endl;
	}
      }
    }
  }

  kseq_destroy(seq);
  gzclose(fp);

  return 0;
}

// Check if all string in the input sample satisfies formulation 3
int check_f3(int argc, char *argv[]) {
  char *index_path = argv[1];
  char *sample_path = argv[2];

  rld_t *index = rld_restore(index_path);

  gzFile fp = gzopen(sample_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  while ((l = kseq_read(seq)) >= 0) {
    uint8_t *P = (uint8_t*)seq->seq.s;
    int m = seq->seq.l;
    seq_char2nt6(seq->seq.l, P);

    bool found_full = backward_search(index, P, m-1);
    bool found_prefix = backward_search(index, P, m-2);
    bool found_suffix = backward_search(index, P+1, m-2);

    if(!found_full && found_prefix && found_suffix)
      ;
    else
      cout << seq->name.s << " " << found_full << found_prefix << found_suffix << endl;
  }
  return 0;
}

vector<vector<pair<string,int>>> analyze_kmers(const vector<pair<string,int>> &kmers, rld_t *index, int threads) {
  vector<vector<pair<string,int>>> solutions (threads);

#pragma omp parallel for num_threads(threads)
  for(uint i=0; i<kmers.size(); ++i) {
    string kmer = kmers[i].first;
    int counter = kmers[i].second;

    char *ckmer = new char[kmer.size()];
    strcpy(ckmer, kmer.c_str());
    uint8_t *cckmer = (uint8_t*)ckmer;
    seq_char2nt6(kmer.size(), cckmer);

    bool hit = backward_search(index, cckmer, kmer.size()-1);

    delete[] ckmer;

    if(!hit)
      solutions[omp_get_thread_num()].push_back(make_pair(kmer,counter));
  }
  return solutions;
}

// Extract child-specific kmers
int kmers(int argc, char *argv[]) {
  char *index_path = argv[1];
  char *kmc_path = argv[2];
  int threads = atoi(argv[3]);

  rld_t *index = rld_restore(index_path);

  CKMCFile kmer_db;
  kmer_db.OpenForListing(kmc_path);
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c, tot_kmers);
  CKmerAPI kmer_obj(klen);

  string kmer;
  int not_found_idx = 1;
  int n_kmers = 0;
  vector<pair<string,int>> kmers;
  string quality (klen, '~');
  while(kmer_db.ReadNextKmer(kmer_obj, counter)) {
    ++n_kmers;
    kmer_obj.to_string(kmer);

    kmers.push_back(make_pair(kmer,counter));

    if(kmers.size() == 50000) {
      vector<vector<pair<string,int>>> solutions = analyze_kmers(kmers, index, threads);

      for(const auto &solution : solutions) {
	for(uint i = 0; i<solution.size(); ++i) {
	  cout << "@kmer_" << not_found_idx << "_" << solution[i].second << endl
	       << solution[i].first << endl
	       << "+" << endl
	       << quality << endl;
	  ++not_found_idx;
	}
      }

      cerr << n_kmers << "/" << tot_kmers << "\r" << flush;
      kmers.clear();
    }
  }

  if(!kmers.empty()) {
    vector<vector<pair<string,int>>> solutions = analyze_kmers(kmers, index, threads);

    for(const auto &solution : solutions) {
      for(uint i = 0; i<solution.size(); ++i) {
	cout << "@kmer_" << not_found_idx << "_" << solution[i].second << endl
	     << solution[i].first << endl
	     << "+" << endl
	     << quality << endl;
	++not_found_idx;
      }
    }
  }
  cerr << tot_kmers << "/" << tot_kmers << endl;

  return 0;
}

int main(int argc, char *argv[]) {
  string mode = argv[1];
  int retcode = 0;
  if(mode == "index")
    retcode = main_build(argc-1, argv+1);
  else if(mode == "sf3")
    retcode = search_f3(argc-1, argv+1);
  else if(mode == "cf3")
    retcode = check_f3(argc-1, argv+1);
  else if(mode == "kmer")
    retcode = kmers(argc-1, argv+1);
  else
    retcode = 1;

  // char *index_path = argv[2];
  // char *read = argv[3];

  // rld_t *index = rld_restore(index_path);

  // uint8_t *P = (uint8_t*)read;
  // int m = strlen(read);
  // int p2 = m - 1;

  // seq_char2nt6(m, P); // convert to integers

  // fmintv_t sai;
  // fm6_set_intv(index, P[p2], sai);
  // cerr << int2char[P[p2]] << " " << interval2str(sai) << endl;
  // --p2;
  // fmintv_t osai[6]; // output SA intervals (one for each symbol between 0 and 5)
  // fm6_extend(index, &sai, osai, 1);
  // for(const auto i : osai)
  //   cerr << interval2str(i) << " ";
  // cerr << endl;
  // sai = osai[P[p2]];
  // cerr << int2char[P[p2]] << " " << interval2str(sai) << endl;

  // fm6_set_intv(index, P[p2], sai);
  // cerr << int2char[P[p2]] << " " << interval2str(sai) << endl;
  // ++p2;
  // fm6_extend(index, &sai, osai, 0);
  // for(const auto i : osai)
  //   cerr << interval2str(i) << " ";
  // cerr << endl;
  // sai = osai[P[p2] >= 1 && P[p2] <= 4 ? 5 - P[p2] : P[p2]];
  // cerr << int2char[P[p2]] << " " << interval2str(sai) << endl;

  return retcode;
}
