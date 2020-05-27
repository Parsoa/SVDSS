#include <iostream>
#include <string>
#include <map>

#include <zlib.h>

#include <omp.h>

#include "kseq.h"
#include "sdsl/suffix_arrays.hpp"
#include "kmc_api/kmc_file.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace sdsl;

typedef csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> FMI;

// TODO: do it better (the current version is too naive)
string reverse_and_complement(const string &S) {
  string rcS;
  map<char,char> encode;
  encode['A'] = 'T';
  encode['C'] = 'G';
  encode['G'] = 'C';
  encode['T'] = 'A';
  for(int i=S.size()-1; i>=0; --i) {
    if(encode.find(S[i]) != encode.end())
      rcS += encode.at(S[i]);
    else {
      rcS += S[i];
      cerr << "Unknown symbol " << S[i] << endl;
    }
  }
  return rcS;
}

// TODO multi-threading
// TODO reduce the number of copy operations
string concatenate_sample(const string &sample_path) {
  gzFile fp = gzopen(sample_path.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  string text = "";
  while ((l = kseq_read(seq)) >= 0)
    text += string(seq->seq.s) + "|" + reverse_and_complement(string(seq->seq.s)) + "|";
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

vector<vector<string>> analyze_kmers(const vector<string> &kmers, const FMI &fm_index, int threads) {
  vector<vector<string>> solutions (threads);

#pragma omp parallel for num_threads(threads)
  for(uint i=0; i<kmers.size(); ++i) {
    string kmer = kmers[i];
    size_t occs = count(fm_index, kmer.begin(), kmer.end());
    if(occs == 0)
      solutions[omp_get_thread_num()].push_back(kmer);
  }
  return solutions;
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
  vector<string> kmers;

  cout << tot_kmers << endl;

  while(kmer_db.ReadNextKmer(kmer_obj, counter)) {
    ++curr_kmer;
    kmer_obj.to_string(kmer);
    kmers.push_back(kmer);

    if(kmers.size() == 10000) {
      vector<vector<string>> unique_kmers = analyze_kmers(kmers, fm_index, threads);

      for(const vector<string> uks : unique_kmers)
	for(const string uk : uks)
	  cout << uk << endl;

      cerr << curr_kmer << "/" << tot_kmers << "\r" << flush;
      kmers.clear();
    }
  }
  if(!kmers.empty()) {
    vector<vector<string>> unique_kmers = analyze_kmers(kmers, fm_index, threads);

    for(const vector<string> uks : unique_kmers)
      for(const string uk : uks)
	cout << uk << endl;
  }
  cerr << curr_kmer << "/" << tot_kmers << endl;

  return 0;
}

int main(int argc, char *argv[]) {
  string mode = argv[1];
  int retcode = 0;
  if(mode == "index")
    retcode = index(argv);
  else if(mode == "count")
    retcode = count(argv);
  else
    retcode = 1;
  return retcode;
}
