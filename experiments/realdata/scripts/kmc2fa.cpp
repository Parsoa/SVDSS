// g++ -std=c++11 -I/miniconda3/include/ -o kmc2fa kmc2fa.cpp -L/miniconda3/lib -lkmc

#include <iostream>

#include "kmc_file.h"

int main(int argc, char *argv[])
{
  char *indbpath = argv[1];
  CKMCFile kmer_db;
  kmer_db.OpenForListing(indbpath);
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c, tot_kmers);
  CKmerAPI kmer_obj(klen);
  char kmer[klen + 1];
  uint i = 0;
  while (kmer_db.ReadNextKmer(kmer_obj, counter)) {
    kmer_obj.to_string(kmer);
    kmer[klen] = '\0';
    std::cout << ">" << "kmer-" << i << "." << counter << "\n"
	      << kmer << std::endl;
    ++i;
    if (i%500 == 0)
      std::cerr << "Parsed " << i << " kmers." << "\r";
  }
  std::cerr << "Parsed " << i << " kmers." << std::endl;
}
