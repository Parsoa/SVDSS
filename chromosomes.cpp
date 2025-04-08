#include "chromosomes.hpp"

KSEQ_INIT(gzFile, gzread)

// FIXME: avoid this
vector<string> chromosomes;
unordered_map<string, char *> chromosome_seqs;

void load_chromosomes(string path) {
  spdlog::info("Loading reference genome from {}..", path);
  gzFile fp = gzopen(path.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  while ((l = kseq_read(seq)) >= 0) {
    string name(seq->name.s);
    spdlog::debug("Extracted {} ({}bp)", seq->name.s, l);
    for (int i = 0; i < l; i++)
      seq->seq.s[i] = toupper(seq->seq.s[i]);
    chromosomes.push_back(seq->name.s);
    char *s = (char *)malloc(sizeof(char) * (l + 1));
    memcpy(s, seq->seq.s, l + 1);
    s[l] = '\0';
    chromosome_seqs[seq->name.s] = s;
  }
  kseq_destroy(seq);
  gzclose(fp);
}

void destroy_chromosomes() {
  for(const auto &chrom : chromosomes) {
    free(chromosome_seqs[chrom]);
  }
}
