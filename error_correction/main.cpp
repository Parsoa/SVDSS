#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <assert.h>

#include <zlib.h>

#include "kseq.h"
#include "kmc_api/kmc_file.h"

#include "bloomfilter.hpp"

KSEQ_INIT(gzFile, gzread)

using namespace std;

auto start_t = chrono::high_resolution_clock::now();

uint8_t char2int[128] = {5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
 			 5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
 			 5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
 			 5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
 			 5, 0, 5, 1,  5, 5, 5, 2,  5, 5, 5, 5,  5, 5, 5, 5,
 			 5, 5, 5, 5,  3, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
 			 5, 0, 5, 1,  5, 5, 5, 2,  5, 5, 5, 5,  5, 5, 5, 5,
 			 5, 5, 5, 5,  3, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5};



void pelapsed(const string &f = "-", const string &s = "-", const bool rollback = false) {
  auto now_t = chrono::high_resolution_clock::now();
  cerr << "[" << f << "/" << s << "] Time elapsed "
	    << chrono::duration_cast<chrono::milliseconds>(now_t - start_t).count()/1000 << "s";
  if(rollback) cerr << "\r";
  else cerr << endl;
}


uint64_t encode(char *s, uint l) {
  // 0<=strlen(s)<=31
  uint64_t n = 0;
  for(uint i=0; i<l; ++i) {
    // i < 128
    if(char2int[s[i]] > 3) return -1;
    n = (n<<2) | char2int[s[i]];
  }
  return n;
}

uint64_t encode_rc(char *s, uint l) {
  uint64_t n = 0;
  for(int i=l-1; i>=0; --i) {
    if(char2int[s[i]] > 3) return -1;
    n = (n<<2) | (3-char2int[s[i]]);
  }
  return n;
}

int main(int argc, char* argv[]) {
  char *fq_path = argv[1];
  char *kmc_prefix = argv[2];
  uint8_t min = atoi(argv[3]);

  pelapsed(__func__, "Initializing BF");
  BF<((uint32_t)0b1 << 31), 3> *bf = new BF<((uint32_t)0b1 << 31), 3>();

  pelapsed(__func__, "Filling BF");
  CKMCFile kmer_db;
  kmer_db.OpenForListing(kmc_prefix);
  uint32 k, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(k, mode, min_counter, pref_len, sign_len, min_c, max_c, tot_kmers);
  CKmerAPI kmer_obj(k);

  char* kmer = new char[k];
  uint64_t enc_kmer = 0;
  uint64_t enc_rckmer = 0;
  while(kmer_db.ReadNextKmer(kmer_obj, counter)) {
    // we assume canonical kmers from KMC (no -b option)
    if(counter < min) continue;
    kmer_obj.to_string(kmer);
    enc_kmer = encode(kmer, k);
    if(enc_kmer != (uint64_t)-1) bf->add(enc_kmer);
  }

  pelapsed(__func__, "Splitting reads");
  gzFile fp = gzopen(fq_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  vector<int> starts;
  vector<int> ends;
  while ((l = kseq_read(seq)) >= 0) {
    for(uint i=0; i<seq->seq.l-k+1; ++i) {
      enc_kmer = encode(seq->seq.s+i, k);
      enc_rckmer = encode_rc(seq->seq.s+i, k);
      enc_kmer = enc_kmer < enc_rckmer ? enc_kmer : enc_rckmer;
      if(starts.size() == ends.size()) {
	// Looking for a start: if the k-mer is not in the BF, we
	// consider the next one, otherwise, we store a start
	if(enc_kmer != (uint64_t)-1 && bf->test(enc_kmer))
	  starts.push_back(i);
      } else {
	// Looking for an end: if the k-mer is in the BF, we consider
	// the next one, otherwise, we store an end and we skip k-1
	// kmers
	if(enc_kmer == (uint64_t)-1 || !bf->test(enc_kmer)) {
	  i+=k-1;
	  ends.push_back(i);
	}
      }
    }
    if(starts.size() != ends.size())
      ends.push_back(seq->seq.l); // end of seq
    assert(starts.size() == ends.size());

    int st,end;
    for(uint i=0; i<starts.size(); ++i) {
      st = starts[i];
      end = ends[i];
      seq->seq.s[end] = '\0';
      cout << ">" << seq->name.s << "_" << st << ":" << end << endl
	   << seq->seq.s+st << endl;
    }
    starts.clear();
    ends.clear();
  }

  delete bf;

  kseq_destroy(seq);
  gzclose(fp);

  pelapsed(__func__, "Done");

  return 0;
}
