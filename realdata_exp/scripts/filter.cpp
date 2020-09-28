// g++ -std=c++11 -O3 filter_specific.cpp -o filter_specific -lz

#include <iostream>

#include <zlib.h>

#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

int filter_by_abundance(int argc, char *argv[]) {
  char* fq_path = argv[1];
  int min_ab = atoi(argv[2]);

  gzFile fp = gzopen(fq_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l, ab;
  while ((l = kseq_read(seq)) >= 0) {
    char *token = strtok(seq->name.s, "#");
    token = strtok(NULL, "#");
    ab = atoi(token);
    if(ab >= min_ab)
      cout << "@" << seq->name.s << "#" << token << "\n"
	   << seq->seq.s << "\n"
	   << "+" << "\n"
	   << seq->qual.s << "\n";
  }
  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}

int filter_by_length(int argc, char *argv[]) {
  char* fq_path = argv[1];
  int max_len = atoi(argv[2]);

  gzFile fp = gzopen(fq_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  while ((l = kseq_read(seq)) >= 0) {
    if(l <= max_len)
      cout << "@" << seq->name.s << "\n"
	   << seq->seq.s << "\n"
	   << "+" << "\n"
	   << seq->qual.s << "\n";
  }
  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}

int filter_by_both(int argc, char *argv[]) {
  char* fq_path = argv[1];
  int min_ab = atoi(argv[2]);
  int max_len = atoi(argv[3]);

  gzFile fp = gzopen(fq_path, "r");
  kseq_t *seq = kseq_init(fp);
  int l, ab;
  while ((l = kseq_read(seq)) >= 0) {
    if(l <= max_len) {
      char *token = strtok(seq->name.s, "#");
      token = strtok(NULL, "#");
      ab = atoi(token);
      if(ab >= min_ab)
	cout << "@" << seq->name.s << "#" << token << "\n"
	     << seq->seq.s << "\n"
	     << "+" << "\n"
	     << seq->qual.s << "\n";
    }
  }
  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}

int main(int argc, char *argv[]) {
  int retcode = 1;
  if(strcmp(argv[1], "l") == 0)
    retcode = filter_by_length(argc-1, argv+1);
  else if(strcmp(argv[1], "a") == 0)
    retcode = filter_by_abundance(argc-1, argv+1);
  else if(strcmp(argv[1], "la") == 0 || strcmp(argv[1], "al") == 0)
    retcode = filter_by_both(argc-1, argv+1);

  return retcode;
}
