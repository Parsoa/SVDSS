#include <iostream>
#include <string>

#include <zlib.h>

#include <omp.h>

#include "kseq.h"
#include "sdsl/suffix_arrays.hpp"
#include "kmc_api/kmc_file.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace sdsl;

typedef csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> FMI;

void build_fmi(const string &sample_path, FMI &fm_index) {
  string index_path = sample_path + ".fm9";
  if (!load_from_file(fm_index, index_path)) {
    cerr << "Building FMI for sample " << sample_path << endl;
    gzFile fp = gzopen(sample_path.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    int l;
    string text = "";
    while ((l = kseq_read(seq)) >= 0)
      text += string(seq->seq.s) + "|";
    kseq_destroy(seq);
    gzclose(fp);
    construct_im(fm_index, text, 1);
    store_to_file(fm_index, index_path);
    cerr << "FMI stored to " << index_path << endl;
  } else {
    cerr << "FMI read from  " << index_path << endl;
  }
}

int analyze_kmers(const vector<string> &kmers, const FMI &fm_index1, const FMI &fm_index2, int threads) {
  vector<int> solutions (threads);

#pragma omp parallel for num_threads(threads)
  for(uint i=0; i<kmers.size(); ++i) {
    string kmer = kmers[i];
    size_t occs1 = count(fm_index1, kmer.begin(), kmer.end());
    size_t occs2 = count(fm_index2, kmer.begin(), kmer.end());
    if(occs1 == 0 && occs2 == 0)
      ++solutions[omp_get_thread_num()];
  }
  int unique_kmers = 0;
  for(const int s : solutions)
    unique_kmers+=s;
  return unique_kmers;
}

int main(int argc, char *argv[]) {
  string sample1_path = argv[1];
  string sample2_path = argv[2];
  string kmc_path = argv[3];
  int threads = stoi(argv[4]);

  FMI fm_index1;
  build_fmi(sample1_path, fm_index1);

  FMI fm_index2;
  build_fmi(sample2_path, fm_index2);

  cerr << "FMI1: " << size_in_mega_bytes(fm_index1) << " MiB" << endl;
  cerr << "FMI2: " << size_in_mega_bytes(fm_index2) << " MiB" << endl;

  CKMCFile kmer_db;
  kmer_db.OpenForListing(kmc_path);
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c, tot_kmers);
  CKmerAPI kmer_obj(klen);
  string kmer;
  int curr_kmer = 0;
  int unique_kmers = 0;
  vector<string> kmers;
  while(kmer_db.ReadNextKmer(kmer_obj, counter)) {
    ++curr_kmer;
    kmer_obj.to_string(kmer);
    kmers.push_back(kmer);

    if(kmers.size() == 10000) {
      int s = analyze_kmers(kmers, fm_index1, fm_index2, threads);
      unique_kmers += s;
      cerr << curr_kmer << "/" << tot_kmers << "\r" << flush;
      kmers.clear();
    }
  }
  if(!kmers.empty()) {
    int s = analyze_kmers(kmers, fm_index1, fm_index2, threads);
    unique_kmers += s;
  }
  cerr << curr_kmer << "/" << tot_kmers << endl;
  cout << unique_kmers << "/" << tot_kmers << endl;

  return 0;
}
