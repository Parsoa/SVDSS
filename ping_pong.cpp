#include "ping_pong.hpp"

// Compute SFS strings from P and store them into solutions
void PingPong::ping_pong_search(const rb3_fmi_t *index, const string &qname,
                                uint8_t *P, int l, vector<SFS> &solutions,
                                int hp_tag) {
  rb3_sai_t ik;
  int begin = l - 1;
  while (begin >= 0) {
    // Backward search. Stop at first mismatch.
    int bmatches = 0;
    rb3_fmd_set_intv(index, P[begin], &ik);

    bmatches = 0;
    while (ik.size != 0 && begin > 0) {
      --begin;
      ++bmatches;
      rb3_sai_t ok[RB3_ASIZE]; // output SA intervals (one for each symbol
                               // between 0 and 5)
      rb3_fmd_extend(index, &ik, ok, 1);
      ik = ok[P[begin]];
    }
    // last checked char (i.e., first of the query) was a match
    if (begin == 0 && ik.size != 0)
      break;

    //  Forward search:
    int end = begin;
    int fmatches = 0;
    rb3_fmd_set_intv(index, P[end], &ik);
    while (ik.size != 0) {
      ++end;
      ++fmatches;
      rb3_sai_t ok[RB3_ASIZE];
      rb3_fmd_extend(index, &ik, ok, 0);
      ik = ok[P[end] >= 1 && P[end] <= 4 ? 5 - P[end] : P[end]];
    }

    int sfs_len = end - begin + 1;
    SFS sfs(qname, begin, sfs_len, hp_tag);
    solutions.push_back(sfs);
    if (begin == 0)
      break;
    if (config->overlap == 0) // Relaxed
      begin -= 1;
    else
      begin = end + config->overlap; // overlap < 0
  }
}

// Load batch from BAM file and store to input entry p. The logic behind is:
// fill position i per each thread, then move to position i+1..
bool PingPong::load_batch_bam(int p) {
  int i = 0;     // current position per thread where to load read
  int nseqs = 0; // loaded seqs

  int o;
  while ((o = sam_read1(bam_file, bam_header,
                        bam_entries[p][nseqs % config->threads][i])) >= 0) {
    bam1_t *alignment = bam_entries[p][nseqs % config->threads][i];
    if (alignment == nullptr) {
      spdlog::critical("nullptr. Why are we here? Please check");
      exit(1);
    }
    reads_processed += 1;
    if (alignment->core.flag & BAM_FUNMAP ||
        alignment->core.flag & BAM_FSUPPLEMENTARY ||
        alignment->core.flag & BAM_FSECONDARY)
      continue;
    if (alignment->core.l_qseq < 100) {
      // FIXME: why do we need this?
      spdlog::warn(
          "Alignment filtered due to l_qseq. Why are we here? Please check");
      continue;
    }
    if (alignment->core.tid < 0) {
      spdlog::critical("core.tid < 0. Why are we here? Please check");
      exit(1);
    }

    uint l = alignment->core.l_qseq; // length of the read
    // if allocated space for read is not enough, extend it
    if (read_seq_max_lengths[p][nseqs % config->threads][i] <
        l) { // FIXME: can't we avoid this? just by allocating *A LOT* per read
      free(read_seqs[p][nseqs % config->threads][i]);
      read_seqs[p][nseqs % config->threads][i] =
          (uint8_t *)malloc(sizeof(char) * (l + 1));
      read_seq_max_lengths[p][nseqs % config->threads][i] = l;
    }
    uint8_t *q = bam_get_seq(alignment);
    for (uint _ = 0; _ < l; _++)
      read_seqs[p][nseqs % config->threads][i][_] =
          seq_nt6_table[(int)seq_nt16_str[bam_seqi(q, _)]];
    read_seqs[p][nseqs % config->threads][i][l] = '\0';

    ++nseqs;
    if (nseqs % config->threads == 0)
      // ith position filled for all threads, move to next one
      ++i;
    if (nseqs == config->batch_size)
      return true;
  }
  assert(o == -1);

  // last batch is incomplete since we reached the end of .bam file
  // we need to fill with nullptr in order to stop processing (we have
  // alignments from previous batches still loaded)

  if (nseqs % config->threads == 0) {
    // we do not have to increment i in this case
  } else {
    // clean remaining threads at position i, then increment i to clear next row
    while (nseqs % config->threads != 0) {
      bam_entries[p][nseqs % config->threads][i] = nullptr;
      ++nseqs;
    }
    ++i;
  }

  // clean next position for all threads
  for (int _ = 0; _ < config->threads; ++_)
    bam_entries[p][_][i] = nullptr;

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
    uint l = fastq_entries[p][n % threads][i].seq.size();
    if (read_seq_max_lengths[p][n % threads][i] < l) {
      free(read_seqs[p][n % threads][i]);
      read_seqs[p][n % threads][i] = (uint8_t *)malloc(sizeof(char) * (l + 1));
      read_seq_max_lengths[p][n % threads][i] = l;
    }
    for (uint _ = 0; _ < l; _++) {
      read_seqs[p][n % threads][i][_] = fastq_entries[p][n % threads][i].seq[_];
    }
    read_seqs[p][n % threads][i][l] = '\0';
    rb3_char2nt6(l, read_seqs[p][n % threads][i]); // convert to integers
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

// Process batch loaded in position p for thread
batch_type_t PingPong::process_batch(const rb3_fmi_t *index, int p,
                                     int thread) {
  batch_type_t solutions;
  // store read id once for all strings to save space, is it worth it?
  if (!bam_mode) {
    for (uint j = 0; j < read_seqs[p][thread].size(); j++) {
      ping_pong_search(
          index, read_names[p][thread][j], read_seqs[p][thread][j],
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
      ping_pong_search(index, qname, read_seqs[p][thread][j], aln->core.l_qseq,
                       solutions[qname], hp_t);
    }
  }
  return solutions;
}

// Output batches until batch b. We keep in memory all batches (empty if already
// output), but some have been already output
void PingPong::output_batch(int b) {
  // FIXME: can we avoid keeping empty batches in memory? maybe unnecessary fix
  for (int i = 0; i < b; i++) {       // for each of the unmerged batches
    for (auto &batch : obatches[i]) { // for each thread in batch
      for (auto &read : batch) {      // for each read in thread
        vector<SFS> assembled_SFSs = read.second;
        if (config->assemble) {
          // Assembler a = Assembler();
          assembled_SFSs = Assembler().assemble(read.second);
        }
        bool is_first = true;
        for (const SFS &sfs : assembled_SFSs) {
          // optimize file output size by not outputting read name for every
          // SFS
          cout << (is_first ? read.first : "*") << "\t" << sfs.qs << "\t"
               << sfs.l << "\t" << sfs.htag << "\t" << endl;
          is_first = false;
        }
      }
      batch.clear();
    }
    obatches[i].clear();
  }
}

// Search for specific strings in input .bam/.fq w.r.t. FMD-Index
int PingPong::search() {
  config = Configuration::getInstance();

  // parse arguments
  spdlog::info("Restoring index..");
  rb3_fmi_t index;
  rb3_fmi_restore(&index, config->index.c_str(), 0);
  if (config->bam != "") {
    bam_file = hts_open(config->bam.c_str(), "r");
    bam_header = sam_hdr_read(bam_file);
    bgzf_mt(bam_file->fp.bgzf, 8, 1);
    bam_mode = 1;
  } else if (config->fastq != "") {
    spdlog::warn("FASTX mode is not optimized (higher running times and larger "
                 "SFSs set).");
    spdlog::warn(
        "If you get a segfault, please consider reducing the batch "
        "size (--bsize) to something <= the number of reads in the sample");
    fastq_file = gzopen(config->fastq.c_str(), "r");
    fastq_iterator = kseq_init(fastq_file);
    bam_mode = 0;
  } else {
    spdlog::critical("No .bam/.fa/.fq file provided, aborting..");
    exit(1);
  }
  // allocate all necessary stuff
  int p = 0; // p can be 0 or 1, used to access the entries/batches currently
             // loaded/analyzed
  for (int i = 0; i < 2; i++) {
    // entries vector contains current and next input
    if (bam_mode) {
      bam_entries.push_back(vector<vector<bam1_t *>>(config->threads));
      for (int j = 0; j < config->threads; ++j)
        for (int k = 0; k < config->batch_size / config->threads;
             ++k) // NOTE: doing so, batch size must be > than #threads
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
        vector<vector<uint>>(config->threads)); // current and next output
    for (int j = 0; j < config->threads; j++) {
      for (int k = 0; k < config->batch_size / config->threads; k++) {
        read_seqs[i][j].push_back((uint8_t *)malloc(sizeof(uint8_t) * (30001)));
        read_names[i][j].push_back("");
        read_seq_max_lengths[i][j].push_back(30000);
      }
    }
  }

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

  if (bam_mode)
    loaded_last_batch = !load_batch_bam(p);
  else
    loaded_last_batch =
        !load_batch_fastq(config->threads, config->batch_size, p);

  uint64_t total_sfs = 0;
  int64_t sfs_to_output = 0;
  spdlog::info("Extracting SFS strings on {} threads..", config->threads);

  while (should_process) {
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
          if (sfs_to_output >= config->max_output) {
            output_batch(obatches.size() - 1);
            sfs_to_output = 0;
          }
        }
      } else {
        // other threads, process the batch
        obatches.back()[t - 2] = process_batch(&index, p, t - 2);
      }
    }

    p = (p + 1) % 2;
    if (should_load)
      obatches.push_back(vector<batch_type_t>(config->threads));
    else
      should_process = false;

    time(&curr_time);
    if (curr_time - start_time == 0)
      ++curr_time;
    time(&curr_time);

    // cerr << "Reads loaded so far: " << reads_processed
    //      << ". Reads processed per second: "
    //      << reads_processed / (curr_time - start_time)
    cerr << "Time: " << curr_time - start_time << "\r";
  }
  cerr << endl;

  output_batch(obatches.size());

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
