#include "aligner.hpp"

KSEQ_INIT(gzFile, gzread)

void Aligner::run() {
  auto c = Configuration::getInstance();
  string fa_path = c->reference;
  string bam_path = c->bam;
  int num_batches = c->aggregate_batches;
  string wd = c->workdir;
  int tau = c->cutoff;

  gzFile fa_in = gzopen(fa_path.c_str(), "r");
  kseq_t *reference = kseq_init(fa_in);
  map<string, string> ref;
  int l;
  while ((l = kseq_read(reference)) >= 0)
    ref[reference->name.s] = reference->seq.s;
  kseq_destroy(reference);
  gzclose(fa_in);
  cout << "Tau is " << tau << endl ;

  lprint({"Aligning high-abundance strings from", to_string(num_batches),
          "batches.."});

  // TODO: do we want to make this step parallel?
  map<string, vector<SFS>> SFSs;
  for (int j = 0; j < num_batches; j++) {
    string s_j = std::to_string(j);
    string inpath = c->workdir + "/solution_batch_" + s_j + ".assembled.sfs";
    map<string, vector<SFS>> batchSFSs = parse_sfsfile(inpath, tau);
    cout << "Loaded " << batchSFSs.size() << " SFS from " << inpath << "." << endl ;
    for (map<string, vector<SFS>>::iterator it = batchSFSs.begin();
         it != batchSFSs.end(); ++it) {
      SFSs[it->first] = it->second;
    }
  }

  // Output SAM file
  string outpath = c->workdir + "/alignments.sam";
  ofstream outf(outpath);

  // 3. Parsing bam
  samFile *bam = hts_open(bam_path.c_str(), "r");
  bam_hdr_t *bamhdr = sam_hdr_read(bam);
  outf << bamhdr->text; // sam_hdr_str(bamhdr) was not declared in this scope
                        // (we need header.o)

  // *** Adding Read Groups to SAM header
  // for (auto x = RGs.begin(); x != RGs.end(); ++x)
  //     cout << "@RG\tID:" << x->second << "\tSM:" << x->first << endl;
  // ***
  bam1_t *aln = bam_init1();
  uint naln = 0;
  while (sam_read1(bam, bamhdr, aln) > 0) {
    char *chr = bamhdr->target_name[aln->core.tid];
    char *qname = bam_get_qname(aln);
    uint qlen = aln->core.l_qseq;
    if (SFSs.find(qname) == SFSs.end())
      continue;

    if (aln->core.flag & BAM_FUNMAP) {
      continue;
    }
    if (aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY)
      continue;

    uint8_t *q = bam_get_seq(aln);
    char *qseq = (char *)malloc(qlen + 1);
    for (uint i = 0; i < qlen; ++i)
      qseq[i] = seq_nt16_str[bam_seqi(q, i)];
    qseq[qlen] = '\0';

    q = bam_get_qual(aln);
    char *qqual = (char *)malloc(qlen + 1);
    for (uint i = 0; i < qlen; ++i)
      // FIXME sometimes qual is empty - didn't know why yet
      qqual[i] = ']'; // char(int(q[i]) + 33);
    qqual[qlen] = '\0';

    bool is_rev = bam_is_rev(aln);
    uint flag = 0;
    if (is_rev)
      flag = 16;

    vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
    // 0-based
    // SNPs: (p,t)
    // Deletions (bases only on ref - contiguous on read): no entry
    // Insertions/Clips: (p, -1)
    for (SFS &sfs : SFSs.at(qname)) {
      if (is_rev && !sfs.isreversed) {
        // CHECKME hard clipping may reduce read lenght and invalidate this.
        // FIX: we require soft-clipping (as pbsv does)
        sfs.reverse(qlen);
      }
      vector<pair<int, int>> local_alpairs;
      for (const pair<int, int> &pos : alpairs)
        if (pos.first != -1 && sfs.s <= pos.first && pos.first < sfs.s + sfs.l)
          local_alpairs.push_back(make_pair(pos.first, pos.second));

      // TODO do we need this if?
      if (local_alpairs.empty()) {
        cerr << "EMPTY LOCAL " << qname << " " << sfs.s << " " << sfs.l << endl;
        continue;
      }

      // FILLING STARTING/ENDING -1s:
      // - if clips, we just add pairs til read end
      // - if insertion, we add pairs til first M we can find
      if (local_alpairs.front().second == -1) {
        bool add = false;
        for (int i = alpairs.size() - 1; i >= 0; --i) {
          if (add)
            local_alpairs.insert(local_alpairs.begin(), alpairs.at(i));
          if (alpairs.at(i).second != -1 && add)
            break;
          if (!add && alpairs.at(i).first == local_alpairs.front().first)
            add = true;
        }
      }
      if (local_alpairs.back().second == -1) {
        bool add = false;
        for (int i = 0; i < alpairs.size(); ++i) {
          if (add)
            local_alpairs.push_back(alpairs.at(i));
          if (alpairs.at(i).second != -1 && add)
            break;
          if (!add && alpairs.at(i).first == local_alpairs.back().first)
            add = true;
        }
      }

      uint qs = local_alpairs.front().first;
      uint qe = local_alpairs.back().first;

      // In some (very rare I hope) cases, an insertion follows a deletions (or
      // viceversa). So we are trying to find the first M - that is a non -1 in
      // the pairs - but that pair has -1 on the read
      if (qs == -1 || qe == -1) {
        cerr << "INS-DEL " << qname << "." << sfs.s << ":" << sfs.l << endl;
        continue;
      }

      // If clips, we have trailing -1 in target positions. We have to find the
      // first placed base
      int ts = local_alpairs.front().second;
      if (local_alpairs.front().first == 0 &&
          local_alpairs.front().second == -1)
        // if we have initial clips, we get the position from original
        // alignments
        ts = aln->core.pos;

      if (local_alpairs.front().second == -1 &&
          local_alpairs.back().second == -1) {
        cerr << "FULL CLIP " << qname << "." << sfs.s << ":" << sfs.l << endl;
        continue;
      }

      CIGAR localcigar = rebuild_cigar(ref[chr], qseq, local_alpairs);
      localcigar.fixclips();
      string localqseq(qseq + qs, qe - qs + 1);
      string localqqual(qqual + qs, qe - qs + 1);

      outf << qname << "." << sfs.s << "-" << sfs.s + sfs.l - 1 << "\t"
           << flag << "\t"
           << chr << "\t"
           << ts + 1 << "\t"
           << to_string(aln->core.qual) << "\t"
           << localcigar.to_str() << "\t"
           << "*" << "\t"
           << 0 << "\t"
           << 0 << "\t"
           << localqseq << "\t"
           << localqqual << "\t"
           << "NM:i:" << localcigar.mismatches
           // << "RG:>:" << RGs[qname]
           << endl;
    }
    free(qseq);
    free(qqual);

    ++naln;
    if (naln % 20000 == 0)
      cerr << naln << " alignments done." << endl;
  }
  outf.close();
  bam_destroy1(aln);
  sam_close(bam);
}

// Partial reimplementation of
// https://github.com/pysam-developers/pysam/blob/6ad0a57ef9c9b05d1492e10228ca7bccb5c7b30e/pysam/libcalignedsegment.pyx#L1867
// TODO run additional tests
vector<pair<int, int>> Aligner::get_aligned_pairs(bam1_t *aln) {
  vector<pair<int, int>> result;
  uint pos = aln->core.pos;
  uint qpos = 0;
  uint32_t *cigar = bam_get_cigar(aln);
  for (uint k = 0; k < aln->core.n_cigar; ++k) {
    uint32_t e = *(cigar + k);
    uint op = bam_cigar_op(e);   // e & 15
    uint l = bam_cigar_oplen(e); // e >> 4
    // char opc = bam_cigar_opchr(op);
    if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF) {
      for (uint i = pos; i < pos + l; ++i) {
        result.push_back(make_pair(qpos, i));
        ++qpos;
      }
      pos += l;
    } else if (op == BAM_CINS or op == BAM_CSOFT_CLIP) {
      for (uint i = pos; i < pos + l; ++i) {
        result.push_back(make_pair(qpos, -1));
        ++qpos;
      }
    } else if (op == BAM_CDEL) {
      for (uint i = pos; i < pos + l; ++i) {
        result.push_back(make_pair(-1, i));
      }
      pos += l;
    } else if (op == BAM_CHARD_CLIP) {
      // advances neither
    } else if (op == BAM_CREF_SKIP) {
      for (uint i = pos; i < pos + l; ++i) {
        result.push_back(make_pair(-1, i));
      }
      pos += l;
    } else if (op == BAM_CPAD) {
      // raise NotImplementedError(
      //     "Padding (BAM_CPAD, 6) is currently not supported. "
      //     "Please implement. Sorry about that.")
    }
  }
  return result;
}

// To build the "local" aligned pairs, we introduce gaps (we do not consider any
// pair where p=-1)
CIGAR Aligner::rebuild_cigar(const string &chr, const string &qseq,
                             const vector<pair<int, int>> &alpairs) {
  CIGAR cigar;
  int last_t = alpairs[0].second;
  for (const pair<int, int> pos : alpairs) {
    int p = pos.first;
    int t = pos.second;
    if (last_t != -1 && t != -1) {
      if (last_t != t && last_t != t - 1) {
        // Deletions
        int d = t - last_t - 1;
        cigar.add(d, 'D', d);
      }
    }
    if (p != -1 && t != -1) {
      // (mis)Match
      if (chr[t] != qseq[p])
        cigar.add(1, 'M', 1);
      else
        cigar.add(1, 'M', 0);
    } else if (t == -1) {
      // Insertion (1bp)
      cigar.add(1, 'I', 1);
    }
    last_t = t;
  }
  return cigar;
}
