#include "converter.hpp"

KSEQ_INIT(gzFile, gzread)

void Converter::run() {
  auto c = Configuration::getInstance();
  string fq_path = c->fastq;
  string sfs_path = c->bed;
  int tau = c->cutoff;

  // 1 Parsing and storing sfss
  map<string, vector<SFS>> SFSs = parse_sfsfile(sfs_path, tau);

  // 2 Parsing reads
  map<string, uint> allsfs;
  gzFile fp = gzopen(fq_path.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  uint idx = 0;
  while ((l = kseq_read(seq)) >= 0) {
    string header = seq->name.s;
    if (SFSs.find(header) == SFSs.end()) {
      // cerr << "Skipping " << header << endl;
      continue;
    }
    for (const SFS &sfs : SFSs.at(header)) {
      uint st = sfs.s;
      uint len = sfs.l;
      uint c = sfs.c;
      string sfsseq(seq->seq.s, st, len);

      string rcsfsseq(sfsseq);
      for (uint j = 0; j < sfsseq.size(); ++j)
        rcsfsseq[j] = RCN[sfsseq[j]];
      reverse(rcsfsseq.begin(), rcsfsseq.end());
      sfsseq = min(sfsseq, rcsfsseq);

      string qual(len, '-');
      cout << "@" << header << "#" << idx << "#" << st << "#" << st + len - 1
           << "#" << c << "\n"
           << sfsseq << "\n"
           << "+"
           << "\n"
           << qual << endl;
      ++idx;
    }
  }

  kseq_destroy(seq);
  gzclose(fp);
}
