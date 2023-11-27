#include "smoother.hpp"

// typedef struct bam1_t {
//     bam1_core_t core; // won't change
//     uint64_t id;      // won't change
//     uint8_t *data;    // has to reallocate
//     int l_data;       // will change
//     uint32_t m_data;  // won't change
//     uint32_t mempolicy:2, :30 /* Reserved */;
// } bam1_t;

// typedef struct bam1_core_t {
//     hts_pos_t pos;      // won't change
//     int32_t tid;        // won't change
//     uint16_t bin;       // just copy
//     uint8_t qual;       // won't change
//     uint8_t l_extranul; // won't change
//     uint16_t flag;      // won't change
//     uint16_t l_qname;   // won't change
//     uint32_t n_cigar;   // will change
//     int32_t l_qseq;     // will change
//     int32_t mtid;       // just copy
//     hts_pos_t mpos;     // just copy
//     hts_pos_t isize;    // may or may ont change?
// } bam1_core_t;

int do_realloc_bam_data(bam1_t *b, size_t desired) {
  uint32_t new_m_data;
  uint8_t *new_data;
  new_m_data = desired;
  kroundup32(new_m_data);
  if (new_m_data < desired) {
    errno = ENOMEM; // Not strictly true but we can't store the size
    return -1;
  }
  new_data = (uint8_t *)realloc(b->data, new_m_data);
  if (!new_data)
    return -1;
  b->data = new_data;
  b->m_data = new_m_data;
  return 0;
}

int realloc_bam_data(bam1_t *b, size_t desired) {
  if (desired <= b->m_data)
    return 0;
  return do_realloc_bam_data(b, desired);
}

void rebuild_bam_entry(bam1_t *alignment, char *seq, uint8_t *qual,
                       vector<pair<uint32_t, uint32_t>> cigar) {
  auto l_aux = bam_get_l_aux(alignment);
  uint8_t *aux = (uint8_t *)malloc(sizeof(uint8_t) * l_aux);
  memcpy(aux, alignment->data + alignment->l_data - l_aux, l_aux);
  // update core
  alignment->core.n_cigar = cigar.size();
  int l = strlen(seq);
  alignment->core.l_qseq = l;
  // rebuild data
  int l_data = alignment->core.l_qname + (4 * alignment->core.n_cigar) +
               ((l + 1) >> 1) + l + l_aux;
  realloc_bam_data(alignment, l_data);
  alignment->l_data = l_data;
  // copy qname
  int offset = alignment->core.l_qname;
  // copy cigar
  uint8_t *cigar_encoded = encode_cigar(cigar);
  memcpy(alignment->data + offset, cigar_encoded, 4 * alignment->core.n_cigar);
  offset += 4 * alignment->core.n_cigar;
  free(cigar_encoded);
  // copy seq data - have to convert seq
  uint8_t *seq_bytes = encode_bam_seq(seq);
  memcpy(alignment->data + offset, seq_bytes, (l + 1) >> 1);
  free(seq_bytes);
  offset += ((l + 1) >> 1);
  // copy quality
  memcpy(alignment->data + offset, qual, l);
  offset += l;
  // copy aux
  memcpy(alignment->data + offset, aux, l_aux);
  free(aux);
}

void Smoother::smooth_read(bam1_t *alignment, char *read_seq, int _i, int _j,
                           int _k) {
  auto cigar_offsets = decode_cigar(alignment);
  int l = 0;
  // try and filter unintenresting reads early on
  bool should_ignore = true;
  for (auto op : cigar_offsets)
    l += op.first;
  if (new_read_seq_max_lengths[_i][_j][_k] < l) {
    free(new_read_seqs[_i][_j][_k]);
    free(new_read_quals[_i][_j][_k]);
    //
    new_read_seqs[_i][_j][_k] = (char *)malloc(sizeof(char) * (l + 1));
    new_read_quals[_i][_j][_k] = (uint8_t *)malloc(sizeof(char) * (l + 1));
    new_read_seq_max_lengths[_i][_j][_k] = l;
  }
  //
  int n = 0;
  size_t m = 0;
  int ref_offset = alignment->core.pos;
  int ins_offset = 0;
  int del_offset = 0;
  int match_offset = 0;
  int soft_clip_offset = 0;
  char *new_seq = new_read_seqs[_i][_j][_k];
  uint8_t *qual = bam_get_qual(alignment);
  uint8_t *new_qual = new_read_quals[_i][_j][_k];
  // Modify current bam1_t* struct
  vector<pair<uint32_t, uint32_t>> new_cigar;
  int m_diff = 0;
  double num_match = 0;
  double num_mismatch = 0;
  char *ref_seq = chromosome_seqs[bam_header->target_name[alignment->core.tid]];
  while (true) {
    if (m == cigar_offsets.size())
      break;
    if (cigar_offsets[m].second == BAM_CMATCH ||
        cigar_offsets[m].second == BAM_CEQUAL ||
        cigar_offsets[m].second == BAM_CDIFF) {
      memcpy(new_seq + n, ref_seq + ref_offset, cigar_offsets[m].first);
      memcpy(new_qual + n, qual + match_offset + ins_offset + soft_clip_offset,
             cigar_offsets[m].first);
      n += cigar_offsets[m].first;
      for (int j = 0; j < cigar_offsets[m].first; j++) {
        num_mismatch +=
            1 ? ref_seq[ref_offset + j] !=
                    read_seq[match_offset + ins_offset + soft_clip_offset + j]
              : 0;
      }
      num_match += cigar_offsets[m].first;
      ref_offset += cigar_offsets[m].first;
      match_offset += cigar_offsets[m].first;
      if (new_cigar.size() >= 1 &&
          new_cigar[new_cigar.size() - 1].second == BAM_CMATCH)
        new_cigar[new_cigar.size() - 1].first +=
            cigar_offsets[m].first + m_diff;
      else
        new_cigar.push_back(
            make_pair(cigar_offsets[m].first + m_diff, BAM_CMATCH));
      m_diff = 0;
    } else if (cigar_offsets[m].second == BAM_CINS) {
      if (cigar_offsets[m].first <= config->min_indel_length) {
        // if a short INDEL then just don't add it to read
      } else {
        // for long INS, this is probably a SV so add it to the read
        should_ignore = false;
        memcpy(new_seq + n,
               read_seq + soft_clip_offset + match_offset + ins_offset,
               cigar_offsets[m].first);
        memcpy(new_qual + n,
               qual + soft_clip_offset + match_offset + ins_offset,
               cigar_offsets[m].first);
        n += cigar_offsets[m].first;
        new_cigar.push_back(cigar_offsets[m]);
      }
      ins_offset += cigar_offsets[m].first;
    } else if (cigar_offsets[m].second == BAM_CDEL) {
      if (cigar_offsets[m].first <= config->min_indel_length) {
        // if a short DEL so let's just fix it
        memcpy(new_seq + n, ref_seq + ref_offset, cigar_offsets[m].first);
        memcpy(new_qual + n,
               qual + soft_clip_offset + match_offset + ins_offset,
               cigar_offsets[m].first);
        n += cigar_offsets[m].first;
        m_diff += cigar_offsets[m].first;
      } else {
        // for long DEL, this is probably a SV so let it be what it was
        should_ignore = false;
        new_cigar.push_back(cigar_offsets[m]);
      }
      del_offset += cigar_offsets[m].first;
      ref_offset += cigar_offsets[m].first;
    } else if (cigar_offsets[m].second == BAM_CSOFT_CLIP) {
      should_ignore = false;
      memcpy(new_seq + n,
             read_seq + soft_clip_offset + match_offset + ins_offset,
             cigar_offsets[m].first);
      memcpy(new_qual + n, qual + soft_clip_offset + match_offset + ins_offset,
             cigar_offsets[m].first);
      n += cigar_offsets[m].first;
      soft_clip_offset += cigar_offsets[m].first;
      new_cigar.push_back(cigar_offsets[m]);
    } else {
      break;
      // if (cigar_offsets[m].second == BAM_CPAD || cigar_offsets[m].second ==
      // BAM_CHARD_CLIP || cigar_offsets[m].second == BAM_CBACK) {
    }
    m += 1;
  }
  new_seq[n] = '\0';
  new_qual[n] = '\0';

  // char *qname = bam_get_qname(alignment);
  // char op = 'M';
  // cerr << qname << " " << strlen(new_seq) << " " << n << endl;
  // for (const auto p : new_cigar) {
  //   if (p.second == BAM_CMATCH || p.second == BAM_CEQUAL ||
  //       p.second == BAM_CDIFF)
  //     op = 'M';
  //   else if (p.second == BAM_CINS)
  //     op = 'I';
  //   else if (p.second == BAM_CDEL)
  //     op = 'D';
  //   else if (p.second == BAM_CSOFT_CLIP)
  //     op = 'S';
  //   cerr << p.first << op;
  // }
  // cerr << endl;

  if (num_mismatch / num_match >= config->al_accuracy) {
    // read is too dirty or is not interesting, just skip it
    bam_aux_update_int(alignment, "XF", 1);
  } else if (should_ignore) {
    bam_aux_update_int(alignment, "XF", 2);
  } else {
    if ((int)strlen(new_seq) != n) {
      // CHECKME: this now should be fixed (read_seqs must not be freed above -
      // only the new ones) cerr << "[W] " << qname << " " << strlen(new_seq) <<
      // "/" << strlen((char*)new_qual) << " " << n << endl;
      spdlog::warn(
          "This shouldn't happen anymore. If you see this warning, please "
          "open an issue at https://github.com/Parsoa/SVDSS/issues");
      bam_aux_update_int(alignment, "XF", 3);
      return;
    }
    rebuild_bam_entry(alignment, new_seq, new_qual, new_cigar);
    bam_aux_update_int(alignment, "XF", 0);
  }
}

void Smoother::process_batch(vector<bam1_t *> bam_entries, int p, int i) {
  bam1_t *alignment;
  for (size_t b = 0; b < bam_entries.size(); b++) {
    alignment = bam_entries[b];
    if (alignment == nullptr)
      // end of the batch
      break;
    smooth_read(alignment, read_seqs[p][i][b], p, i, b);
  }
}

// BAM writing based on https://www.biostars.org/p/181580/
void Smoother::run() {
  config = Configuration::getInstance();

  load_chromosomes(config->reference);

  // parse arguments
  bam_file = hts_open(config->bam.c_str(), "r");
  bam_index = sam_index_load(bam_file, config->bam.c_str());
  bam_header = sam_hdr_read(bam_file); // read header
  bgzf_mt(bam_file->fp.bgzf, 8, 1);
  out_bam_file = hts_open("-", "wb");
  bgzf_mt(out_bam_file->fp.bgzf, 8, 1);
  if (sam_hdr_write(out_bam_file, bam_header) < 0) {
    spdlog::critical("Can't write corrected BAM header, aborting..");
    return;
  }
  // allocate stuff
  int modulo = 3;
  for (int i = 0; i < modulo; i++) {
    // original read sequences
    read_seqs.push_back(vector<vector<char *>>(config->threads));
    // smoothed read sequences
    new_read_seqs.push_back(vector<vector<char *>>(config->threads));
    // smoothed read qualities
    new_read_quals.push_back(vector<vector<uint8_t *>>(config->threads));
    // original read max lengths
    read_seq_max_lengths.push_back(vector<vector<int>>(config->threads));
    // smoothed read max lengths
    new_read_seq_max_lengths.push_back(vector<vector<int>>(config->threads));
    for (int j = 0; j < config->threads; j++) {
      for (int k = 0; k < config->batch_size / config->threads; k++) {
        read_seqs[i][j].push_back((char *)malloc(sizeof(char) * (30001)));
        new_read_seqs[i][j].push_back((char *)malloc(sizeof(char) * (30001)));
        new_read_quals[i][j].push_back(
            (uint8_t *)malloc(sizeof(uint8_t) * (30001)));
        //
        read_seq_max_lengths[i][j].push_back(30000);
        new_read_seq_max_lengths[i][j].push_back(30000);
      }
    }
  }
  for (int i = 0; i < modulo; i++) {
    bam_entries.push_back(vector<vector<bam1_t *>>(config->threads));
    for (int j = 0; j < config->threads; j++)
      for (int k = 0; k < config->batch_size / config->threads; k++)
        bam_entries[i][j].push_back(bam_init1());
  }

  spdlog::info("Smoothing alignments on {} threads..", config->threads);
  time_t start_time;
  time_t curr_time;
  time(&start_time);
  int b = 0;
  int p = 1;
  load_batch_bam(p);
  bool should_load = true;
  bool should_process = true;
  bool should_terminate = false;
  bool loaded_last_batch = false;
  int reads_written = 0;
  while (should_process) {
    if (!should_load)
      should_process = false;
    if (loaded_last_batch)
      should_load = false;
#pragma omp parallel for num_threads(config->threads + 2)
    for (int i = 0; i < config->threads + 2; i++) {
      if (i == 0) {
        // first thread loads next batch
        if (should_load)
          loaded_last_batch = !load_batch_bam((p + 1) % modulo);
      } else if (i == 1) {
        // second thread output batches
        if (b >= 1) {
          int ret = 0;
          for (int k = 0; k < config->batch_size / config->threads; k++) {
            for (int j = 0; j < config->threads; j++) {
              if (bam_entries[(p + 2) % modulo][j][k] != nullptr) {
                if (bam_entries[(p + 2) % modulo][j][k]->core.flag &
                        BAM_FUNMAP ||
                    bam_entries[(p + 2) % modulo][j][k]->core.flag &
                        BAM_FSUPPLEMENTARY ||
                    bam_entries[(p + 2) % modulo][j][k]->core.flag &
                        BAM_FSECONDARY) {
                  // CHECKME
                  spdlog::warn("We have a non primary smoothed alignment. "
                               "How can we be here?");
                  continue;
                }
                ret = sam_write1(out_bam_file, bam_header,
                                 bam_entries[(p + 2) % modulo][j][k]);
                reads_written += 1;
                if (ret < 0) {
                  spdlog::critical(
                      "Can't write corrected BAM record, aborting..");
                  should_terminate = true;
                  break;
                }
              } else
                break;
            }
          }
        }
      } else {
        // other threads, process the batch
        if (should_process)
          process_batch(bam_entries[p][i - 2], p, i - 2);
      }
    }
    if (should_terminate)
      return;
    p += 1;
    p %= modulo;
    b += 1;

    time(&curr_time);
    if (curr_time - start_time == 0)
      ++curr_time;
    cerr << "Alignments processed so far: " << reads_processed
         << ". Alignments processed per second: "
         << reads_processed / (curr_time - start_time)
         // << ". Alignments wrote: " << reads_written
         << ". Time: " << curr_time - start_time << "\r";
  }
  sam_close(bam_file);
  sam_close(out_bam_file);

  cerr << endl;
}

/* Load batch from BAM file and store to input entry p. The logic behind is:
 * fill position i per each thread, then move to position i+1.. */
bool Smoother::load_batch_bam(int p) {
  int i = 0;     // current position per thread where to load read
  int nseqs = 0; // loaded seqs
  while (sam_read1(bam_file, bam_header,
                   bam_entries[p][nseqs % config->threads][i]) >= 0) {
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
    if (alignment->core.qual < config->min_mapq) {
      continue;
    }
    if (alignment->core.l_qseq < 2) {
      // FIXME: why do we need this?
      spdlog::warn(
          "Alignment filtered due to l_qseq. Why are we here? Please check");
      continue;
    }
    if (alignment->core.tid < 0) {
      spdlog::critical("core.tid < 0. Why are we here? Please check");
      exit(1);
    }
    if (chromosome_seqs.find(bam_header->target_name[alignment->core.tid]) ==
        chromosome_seqs.end()) {
      if (warned_chromosomes.find(
              bam_header->target_name[alignment->core.tid]) ==
          warned_chromosomes.end())
        spdlog::warn(
            "Skipping alignment(s) on {} since it is not present in the "
            "reference",
            bam_header->target_name[alignment->core.tid]);
      warned_chromosomes.insert(bam_header->target_name[alignment->core.tid]);
      continue;
    }

    int l = alignment->core.l_qseq; // length of the read
    // if allocated space for read is not enough, reallocate more
    if (read_seq_max_lengths[p][nseqs % config->threads][i] <
        l) { // FIXME: can we avoid this just by allocating *A LOT* per read?
      free(read_seqs[p][nseqs % config->threads][i]);
      read_seqs[p][nseqs % config->threads][i] =
          (char *)malloc(sizeof(char) * (l + 1));
      read_seq_max_lengths[p][nseqs % config->threads][i] = l;
    }
    uint8_t *q = bam_get_seq(alignment);
    for (int _ = 0; _ < l; _++)
      read_seqs[p][nseqs % config->threads][i][_] =
          seq_nt16_str[bam_seqi(q, _)];
    read_seqs[p][nseqs % config->threads][i][l] = '\0';
    ++nseqs;
    if (nseqs % config->threads == 0)
      ++i;
    if (nseqs == config->batch_size)
      return true;
  }
  // last batch is incomplete since we reached the end of .bam file
  // TODO: can we do like in ping_pong?
  if (nseqs != config->batch_size) {
    for (int j = nseqs % config->threads; j < config->threads; j++)
      for (int _ = i; _ < config->batch_size / config->threads; _++)
        bam_entries[p][j][_] = nullptr;
    for (int j = 0; j < nseqs % config->threads; j++)
      for (int _ = i + 1; _ < config->batch_size / config->threads; _++)
        bam_entries[p][j][_] = nullptr;
  }
  return false;
}
