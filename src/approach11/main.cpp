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

// // Very naive way to convert a string to a binary string
// // TODO: do better
// string encode_3bits(const string &text) {
//   string text_3b (3*text.size(), '0');

//   map<char, string> encode;
//   encode['A'] = "000";
//   encode['C'] = "001";
//   encode['G'] = "010";
//   encode['T'] = "011";

//   int i = 0;
//   for(const char c : text) {
//     if(encode.find(c) == encode.end())
//       text_3b.replace(3*i, 3, "111");
//     else
//       text_3b.replace(3*i, 3, encode.at(c));
//     ++i;
//   }

//   return text_3b;
// }

int analyze_kmers(const vector<string> &kmers, const FMI &fm_index, int threads) {
  vector<int> solutions (threads);

#pragma omp parallel for num_threads(threads)
  for(uint i=0; i<kmers.size(); ++i) {
    string kmer = kmers[i];
    size_t occs = count(fm_index, kmer.begin(), kmer.end());
    if(occs == 0)
      ++solutions[omp_get_thread_num()];
  }
  int unique_kmers = 0;
  for(const int s : solutions)
    unique_kmers+=s;
  return unique_kmers;
}

/**
 * TODO:
 *      - multi-threading
 *      - reduce the number of copy operations (maybe using cstrings)
 **/
string concatenate_sample(const string &sample_path) {
  gzFile fp = gzopen(sample_path.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  string text = "";
  while ((l = kseq_read(seq)) >= 0)
    text += string(seq->seq.s) + "|";
  kseq_destroy(seq);
  gzclose(fp);
  return text;
}

int index(char *argv[]) {
  string sample1_path = argv[2];
  string sample2_path = argv[3];
  string index_path = argv[4];

  FMI fm_index;
  if (!load_from_file(fm_index, index_path)) {
    cerr << "Building FMI for samples " << sample1_path << "," << sample2_path << endl;
    string text = concatenate_sample(sample1_path) + concatenate_sample(sample2_path);
    construct_im(fm_index, text, 1);
    store_to_file(fm_index, index_path);
    cerr << "FMI stored to " << index_path << " (" << size_in_mega_bytes(fm_index) << "MiB)" << endl;
  } else {
    cerr << "FMI found (" << index_path << ")" << endl;
  }

  return 0;
}

int count(char *argv[]) {
  string index_path = argv[2];
  string kmc_path = argv[3];
  int threads = stoi(argv[4]);

  FMI fm_index;
  if (!load_from_file(fm_index, index_path)) {
    cerr << "FMI not found (" << index_path << ")" << endl;
    return 1;
  }
  cerr << "FMI: " << size_in_mega_bytes(fm_index) << " MiB" << endl;

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
      int s = analyze_kmers(kmers, fm_index, threads);
      unique_kmers += s;
      cerr << curr_kmer << "/" << tot_kmers << "\r" << flush;
      kmers.clear();
    }
  }
  if(!kmers.empty()) {
    int s = analyze_kmers(kmers, fm_index, threads);
    unique_kmers += s;
  }
  cerr << curr_kmer << "/" << tot_kmers << endl;
  cout << unique_kmers << "/" << tot_kmers << endl;

  return 0;
}

int main(int argc, char *argv[]) {
  string mode = argv[1];
  int retcode = 0;
  if(mode == "index")
    retcode = index(argv);
  else if(mode == "count")
    retcode = count(argv);

  return retcode;
}
