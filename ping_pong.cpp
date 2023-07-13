#include "ping_pong.hpp"

// string interval2str(rldintv_t sai) {
//   return "[" + to_string(sai.x[0]) + "," + to_string(sai.x[1]) + "," +
//           to_string(sai.x[2]) + "]";
// }

void seq_char2nt6(int l, unsigned char *s) {
  int i;
  for (i = 0; i < l; ++i) {
    s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
  }
}

/* Compute SFS strings from P and store them into solutions*/
void PingPong::ping_pong_search(rld_t *index, uint8_t *P, int l,
                                vector<sfs_type_t> &solutions, int hp_tag) {
  rldintv_t sai;
  int begin = l - 1;
  while (begin >= 0) {
    // Backward search. Stop at first mismatch.
    int bmatches = 0;
    fm6_set_intv(index, P[begin], sai);
    // cerr << "BS from " << int2char[P[begin]] << " (" << begin << "): " <<
    // interval2str(sai) << endl;
    bmatches = 0;
    while (sai.x[2] != 0 && begin > 0) {
      begin--;
      bmatches++;
      rldintv_t
          osai[6]; // output SA intervals (one for each symbol between 0 and 5)
      rld_extend(index, &sai, osai, 1);
      sai = osai[P[begin]];
    }
    // last checked char (i.e., first of the query) was a match
    if (begin == 0 && sai.x[2] != 0) {
      break;
    }
    // cerr << "Mismatch " << int2char[P[begin]] << " (" <<  begin << ").
    // bmatches: " << to_string(bmatches) << endl;

    //  Forward search:
    int end = begin;
    int fmatches = 0;
    fm6_set_intv(index, P[end], sai);
    // cerr << "FS from " << int2char[P[end]] << " (" << end << "): " <<
    // interval2str(sai) << endl;
    while (sai.x[2] != 0) {
      end++;
      fmatches++;
      rldintv_t osai[6];
      rld_extend(index, &sai, osai, 0);
      sai = osai[P[end] >= 1 && P[end] <= 4 ? 5 - P[end] : P[end]];
      // DEBUG(cerr << "- FE with " << int2char[P[end]] << " (" <<  end << "): "
      // << interval2str(sai) << endl ;)
    }
    // DEBUG(cerr << "Mismatch " << int2char[P[end]] << " (" << end << ").
    // fmatches: " << fmatches << endl ;) DEBUG(cerr << "Adding [" << begin <<
    // ", " << end << "]." << endl ;)
    int sfs_len = end - begin + 1;
    SFS sfs = SFS{begin, sfs_len, hp_tag};
    solutions.push_back(sfs);
    if (begin == 0) {
      break;
    }
    if (config->overlap == 0) { // Relaxed
      begin -= 1;
    } else {
      begin = end + config->overlap; // overlap < 0
    }
  }
}

/* Load batch from BAM file and store to input entry p. The logic behind is:
 * fill position i per each thread, then move to position i+1.. */
bool PingPong::load_batch_bam(int p) {
  int threads = config->threads;
  int batch_size = config->batch_size;
  int i = 0; // current position per thread where to load read
  int nseqs = 0;

  int o;
  while ((o = sam_read1(bam_file, bam_header,
                        bam_entries[p][nseqs % threads][i])) >= 0) {
    bam1_t *alignment = bam_entries[p][nseqs % threads][i];
    if (alignment == nullptr) {
      spdlog::critical("nullptr. Why are we here? Please check");
      exit(1);
    }
    reads_processed += 1;
    if (alignment->core.flag & BAM_FUNMAP ||
        alignment->core.flag & BAM_FSUPPLEMENTARY ||
        alignment->core.flag & BAM_FSECONDARY)
      continue;
    if (alignment->core.l_qseq < 100) // FIXME: why do we need this?
      continue;
    if (alignment->core.tid < 0) {
      spdlog::critical("core.tid < 0. Why are we here? Please check");
      exit(1);
    }

    uint l = alignment->core.l_qseq; // length of the read
    // if allocated space for read is not enough, extend it
    if (read_seq_max_lengths[p][nseqs % threads][i] <
        l) { // FIXME: can't we avoid this? just by allocating *A LOT* per read
      free(read_seqs[p][nseqs % threads][i]);
      read_seqs[p][nseqs % threads][i] =
          (uint8_t *)malloc(sizeof(char) * (l + 1));
      read_seq_max_lengths[p][nseqs % threads][i] = l;
    }
    uint8_t *q = bam_get_seq(alignment);
    for (uint _ = 0; _ < l; _++)
      read_seqs[p][nseqs % threads][i][_] =
          seq_nt16_str[bam_seqi(q, _)] < 128
              ? seq_nt6_table[seq_nt16_str[bam_seqi(q, _)]]
              : 5;
    read_seqs[p][nseqs % threads][i][l] = '\0';

    ++nseqs;
    if (nseqs % threads == 0)
      // ith position filled for all threads, move to next one
      ++i;
    if (nseqs == batch_size)
      return true;
  }
  assert(o == -1);

  // we need to fill with nullptr in order to stop processing (we have
  // alignments from previous batches still loaded)

  // clean remaining threads at position i
  while (nseqs % threads != 0) {
    bam_entries[p][nseqs % threads][i] = nullptr;
    ++nseqs;
  }
  ++i;
  // clean next position for all threads
  for (int _ = 0; _ < threads; ++_)
    bam_entries[p][_ % threads][i] = nullptr;

  return false;
}

bool PingPong::load_batch_fastq(int threads, int batch_size, int p) {
  for (int i = 0; i < threads; i++) {
    fastq_entries[p][i].clear();
  }
  int k = 0;
  int i = 0;
  int n = 0;
  while ((k = kseq_read(fastq_iterator)) >= 0) {
    if (fastq_iterator->qual.l) {
      fastq_entries[p][n % threads].push_back(
          fastq_entry_t(fastq_iterator->name.s, fastq_iterator->seq.s,
                        fastq_iterator->qual.s));
    } else {
      fastq_entries[p][n % threads].push_back(
          fastq_entry_t(fastq_iterator->name.s, fastq_iterator->seq.s,
                        fastq_iterator->name.s));
    }
    int l = fastq_entries[p][n % threads][i].seq.size();
    if (read_seq_max_lengths[p][n % threads][i] < l) {
      free(read_seqs[p][n % threads][i]);
      read_seqs[p][n % threads][i] = (uint8_t *)malloc(sizeof(char) * (l + 1));
      read_seq_max_lengths[p][n % threads][i] = l;
    }
    for (int _ = 0; _ < l; _++) {
      read_seqs[p][n % threads][i][_] = fastq_entries[p][n % threads][i].seq[_];
    }
    read_seqs[p][n % threads][i][l] = '\0';
    seq_char2nt6(l, read_seqs[p][n % threads][i]); // convert to integers
    read_names[p][n % threads][i] = fastq_iterator->name.s;
    n += 1;
    if (n % threads == 0) {
      i += 1;
    }
    if (n == batch_size) {
      return true;
    }
  }
  return n != 0 ? true : false;
}

/* Process batch loaded in position p for thread */
batch_type_t PingPong::process_batch(rld_t *index, int p, int thread) {
  batch_type_t solutions;
  // store read id once for all strings to save space, is it worth it?
  if (!bam_mode) {
    for (uint j = 0; j < read_seqs[p][thread].size(); j++) {
      ping_pong_search(
          index, read_seqs[p][thread][j],
          strlen(
              (char *)read_seqs[p][thread]
                               [j]), // FIXME: this may be inefficient. We were
                                     // storing lengths in a vector as sequences
          solutions[read_names[p][thread][j]], 0);
    }
  } else {
    for (uint j = 0; j < bam_entries[p][thread].size(); j++) {
      bam1_t *aln = bam_entries[p][thread][j];
      if (aln == nullptr)
        break;
      char *qname = bam_get_qname(aln);
      int xf_t = bam_aux_get(aln, "XF") != NULL
                     ? bam_aux2i(bam_aux_get(aln, "XF"))
                     : 0;
      int hp_t = bam_aux_get(aln, "HP") != NULL
                     ? bam_aux2i(bam_aux_get(aln, "HP"))
                     : 0;
      if (config->putative and xf_t != 0)
        continue;
      ping_pong_search(index, read_seqs[p][thread][j], aln->core.l_qseq,
                       solutions[qname], hp_t);
    }
  }
  return solutions;
}

/* Output batches until batch b. We keep in memory all batches (empty if already
 * output), but some have been already output*/
void PingPong::output_batch(int b) {
  auto c = Configuration::getInstance();
  string path = c->workdir + "/solution_batch_" + to_string(current_obatch) +
                (!c->assemble ? ".unassembled" : "") + ".sfs";
  ofstream of(path);
  // FIXME: can we avoid keeping empty batches in memory? maybe unnecessary fix
  for (int i = 0; i < b; i++) {       // for each of the unmerged batches
    for (auto &batch : obatches[i]) { // for each thread in batch
      for (auto &read : batch) {      // for each read in thread
        vector<SFS> assembled_SFSs = read.second;
        if (c->assemble) {
          // Assembler a = Assembler();
          assembled_SFSs = Assembler().assemble(read.second);
        }
        bool is_first = true;
        for (const SFS &sfs : assembled_SFSs) {
          // optimize file output size by not outputting read name for every
          // SFS
          of << (is_first ? read.first : "*") << "\t" << sfs.s << "\t" << sfs.l
             << "\t" << sfs.htag << "\t" << endl;
          is_first = false;
        }
      }
      batch.clear();
    }
    obatches[i].clear();
  }
}

/** Search for specific strings in input .bam/.fq w.r.t. FMD-Index **/
int PingPong::search() {
  config = Configuration::getInstance();

  // parse arguments
  spdlog::info("Restoring index..");
  rld_t *index = rld_restore(config->index.c_str());
  if (config->bam != "") {
    bam_file = hts_open(config->bam.c_str(), "r");
    bam_header = sam_hdr_read(bam_file);
    bgzf_mt(bam_file->fp.bgzf, 8, 1);
    bam_mode = 1;
  } else if (config->fastq != "") {
    spdlog::warn("FASTQ mode is not optimized");
    fastq_file = gzopen(config->fastq.c_str(), "r");
    fastq_iterator = kseq_init(fastq_file);
    bam_mode = 0;
  } else {
    spdlog::critical("No .bam/.fq file provided, aborting..");
    exit(1);
  }
  // allocate all necessary stuff
  int p = 0; // p can be 0 or 1, used to access the entries/batches currently
             // loaded/analyzed
  for (int i = 0; i < 2; i++) {
    // entries vector contains current and next input
    if (bam_mode) {
      bam_entries.push_back(vector<vector<bam1_t *>>(config->threads));
      for (int j = 0; j < config->threads; j++)
        for (int k = 0; k < config->batch_size / config->threads; k++)
          bam_entries[i][j].push_back(bam_init1());
    } else
      fastq_entries.push_back(vector<vector<fastq_entry_t>>(config->threads));
  }
  // pre-allocate read seqs
  for (int i = 0; i < 2; i++) {
    read_seqs.push_back(
        vector<vector<uint8_t *>>(config->threads)); // current and next output
    read_names.push_back(
        vector<vector<string>>(config->threads)); // current and next output
    read_seq_max_lengths.push_back(
        vector<vector<int>>(config->threads)); // current and next output
    for (int j = 0; j < config->threads; j++) {
      for (int k = 0; k < config->batch_size / config->threads; k++) {
        read_seqs[i][j].push_back((uint8_t *)malloc(sizeof(uint8_t) * (30001)));
        read_names[i][j].push_back("");
        read_seq_max_lengths[i][j].push_back(30000);
      }
    }
  }

  if (bam_mode)
    load_batch_bam(p);
  else
    load_batch_fastq(config->threads, config->batch_size, p);

  obatches.push_back(vector<batch_type_t>(
      config->threads)); // each loaded entry will produce a batch. An output
                         // batch is a vector of config->threads batches. We
                         // keep all output batches in memory (empty if already
                         // output) but only two input entries

  // main loop
  time_t start_time;
  time_t curr_time;
  time(&start_time);

  bool should_load = true;
  bool should_process = true;
  bool loaded_last_batch = false;

  uint64_t total_sfs = 0;
  int64_t sfs_to_output = 0;
  spdlog::info("Extracting SFS strings on {} threads..", config->threads);

  while (should_process) {
    if (!should_load)
      should_process = false;
    if (loaded_last_batch)
      should_load = false;
#pragma omp parallel for num_threads(config->threads + 2) schedule(static, 1)
    for (int i = 0; i < config->threads + 2; i++) {
      int t = omp_get_thread_num();
      if (t == 0) {
        // first thread loads next batch
        if (should_load) {
          if (bam_mode)
            loaded_last_batch = !load_batch_bam((p + 1) % 2);
          else
            loaded_last_batch = !load_batch_fastq(
                config->threads, config->batch_size, (p + 1) % 2);
        }
      } else if (t == 1) {
        // second thread counts and output batches (from last one not output yet
        // to last finished)
        if (obatches.size() > 1) {
          uint64_t c = 0;
          for (const auto &batch : obatches[obatches.size() - 2])
            for (auto it = batch.begin(); it != batch.end(); it++)
              c += it->second.size();
          total_sfs += c;
          sfs_to_output += c;
          if (sfs_to_output >= config->max_output || !should_process) {
            if (!should_process)
              output_batch(obatches.size());
            else
              output_batch(obatches.size() - 1);
            sfs_to_output = 0;
            current_obatch += 1;
          }
        }
      } else {
        // other threads, process the batch
        if (should_process)
          obatches.back()[t - 2] = process_batch(index, p, t - 2);
      }
    }

    p = (p + 1) % 2;
    if (should_load)
      obatches.push_back(vector<batch_type_t>(config->threads));

    if (curr_time - start_time == 0)
      ++curr_time;
    time(&curr_time);

    cerr << "Reads loaded so far: " << reads_processed
         << ". Reads processed per second: "
         << reads_processed / (curr_time - start_time)
         << ". Batches exported: " << current_obatch
         << ". Time: " << curr_time - start_time << "\r";
  }
  cerr << endl;

  // cleanup
  if (bam_mode) {
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < config->threads; j++)
        for (int k = 0; k < config->batch_size / config->threads; k++)
          bam_destroy1(bam_entries[i][j][k]);
  } else {
    kseq_destroy(fastq_iterator);
    gzclose(fastq_file);
  }

  return 0;
}

/* Build FMD-index for input .fa/.fq. Code adapted from ropebwt2 (main_ropebwt2
 * in main.c) **/
int PingPong::index() {
  // hardcoded parameters
  uint64_t m = (uint64_t)(.97 * 10 * 1024 * 1024 * 1024) +
               1; // batch size for multi-string indexing
  int block_len = ROPE_DEF_BLOCK_LEN, max_nodes = ROPE_DEF_MAX_NODES,
      so = MR_SO_RCLO;
  int thr_min =
      4; // switch to single thread when < 100 strings remain in a batch

  // the index
  mrope_t *mr = 0;

  bool binary_output = config->binary;
  if (config->append != "") {
    FILE *fp;
    spdlog::info("Appending to index: {}", config->append);
    if ((fp = fopen(config->append.c_str(), "rb")) == 0) {
      spdlog::critical("Failed to open file {}", config->append);
      return 1;
    }
    mr = mr_restore(fp);
    fclose(fp);
  }

  // Initialize mr if not restored
  if (mr == 0)
    mr = mr_init(max_nodes, block_len, so);
  mr_thr_min(mr, thr_min);

  // Parsing the input sample
  gzFile fp;
  if (config->reference != "")
    fp = gzopen(config->reference.c_str(), "rb");
  else if (config->fastq != "")
    fp = gzopen(config->fastq.c_str(), "rb");
  else {
    spdlog::critical("Please provide a FASTA/Q file. Halting..");
    exit(1);
  }
  kseq_t *ks = kseq_init(fp);
  kstring_t buf = {0, 0, 0}; // buffer, will contain the concatenation
  int l;
  uint8_t *s;
  int i;
  while ((l = kseq_read(ks)) >= 0) {
    s = (uint8_t *)ks->seq.s;

    // change encoding
    for (i = 0; i < l; ++i) {
      s[i] = s[i] < 128 ? seq_nt6_table[s[i]] : 5;
    }

    // Reverse the sequence
    for (i = 0; i < (l >> 1); ++i) {
      int tmp = s[l - 1 - i];
      s[l - 1 - i] = s[i];
      s[i] = tmp;
    }

    // Add forward to buffer
    kputsn((char *)ks->seq.s, ks->seq.l + 1, &buf);

    // Add reverse to buffer
    for (i = 0; i < (l >> 1); ++i) {
      int tmp = s[l - 1 - i];
      tmp = (tmp >= 1 && tmp <= 4) ? 5 - tmp : tmp;
      s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
      s[i] = tmp;
    }
    if (l & 1)
      s[i] = (s[i] >= 1 && s[i] <= 4) ? 5 - s[i] : s[i];
    kputsn((char *)ks->seq.s, ks->seq.l + 1, &buf);

    if (buf.l >= m) {
      mr_insert_multi(mr, buf.l, (uint8_t *)buf.s, 1);
      buf.l = 0;
    }
  }

  if (buf.l) { // last batch
    mr_insert_multi(mr, buf.l, (uint8_t *)buf.s, 1);
  }

  free(buf.s);
  kseq_destroy(ks);
  gzclose(fp);

  // dump index to stdout
  if (binary_output) {
    // binary FMR format
    mr_dump(mr, fopen(config->index.c_str(), "wb"));
  } else {
    // FMD format
    mritr_t itr;
    const uint8_t *block;
    rld_t *e = 0;
    rlditr_t di;
    e = rld_init(6, 3);
    rld_itr_init(e, &di, 0);
    mr_itr_first(mr, &itr, 1);
    while ((block = mr_itr_next_block(&itr)) != 0) {
      const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
      while (q < end) {
        int c = 0;
        int64_t l;
        rle_dec1(q, c, l);
        rld_enc(e, &di, l, c);
      }
    }
    rld_enc_finish(e, &di);
    rld_dump(e, config->index.c_str());
  }

  mr_destroy(mr);

  return 0;
}
