#include "realigner.hpp"
#include "chromosomes.hpp"

bool Realigner::load_batch_bam(int threads, int batch_size, int p) {
    int i = 0 ;
    int n = 0 ;
    while (sam_read1(bam_file, bam_header, bam_entries[p][n % threads][i]) >= 0) {
        auto alignment = bam_entries[p][n % threads][i] ;
        if (alignment == nullptr) {
            break ;
        }
        // these reads have not been reconstructed
        if (alignment->core.flag & BAM_FUNMAP || alignment->core.flag & BAM_FSUPPLEMENTARY || alignment->core.flag & BAM_FSECONDARY) {
            continue ;
        }
        if (alignment->core.l_qseq <= 10) {
            //cerr << "Read too short, ignoring.." << endl ;
            continue ;
        }
        if (alignment->core.tid < 0) {
            continue ;
        }
        n += 1 ;
        if (n % threads == 0) {
            i += 1 ;
        }
        if (n == batch_size) {
            break ;
        }
    }
    lprint({"Loaded", to_string(n), "BAM reads.."});
    return n != 0 ? true : false ;
}

// load SFS file and round-robin reads between threads
void Realigner::load_input_sfs_batch() {
    if (current_input_batch == config->aggregate_batches) {
        cerr << "[I] all input SFS batches loaded. Not loading new batch." << endl ;
        return ;
    }
    string s_j = std::to_string(current_input_batch) ;
    string inpath = config->workdir + "/solution_batch_" + s_j + ".assembled.sfs" ;
    cout << "[I] Loading SFS from " << inpath << endl ;
    int threads = config->threads ; 
    unordered_map<string, vector<SFS>> SFSs ;
    string line ;
    ifstream inf(inpath) ;
    if (inf.is_open()) {
        string info[4];
        string read_name;
        while (getline(inf, line)) {
            stringstream ssin(line);
            int i = 0;
            while (ssin.good() && i < 5) {
                ssin >> info[i++];
            }
            if (info[0].compare("*") != 0) {
                read_name = info[0];
                SFSs[read_name] = vector<SFS>();
            }
            SFSs[read_name].push_back(SFS(stoi(info[1]), stoi(info[2]), stoi(info[3]), true)) ;
        }
    }
    sfs_batches.push_back(SFSs) ;
    current_input_batch++ ;
}

void Realigner::run() {
    config = Configuration::getInstance() ;
    load_chromosomes(config->reference) ;
    bam_file = sam_open(config->bam.c_str(), "r") ;
    bam_header = sam_hdr_read(bam_file) ; //read header
    // load first batch
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(config->threads)) ;
    }
    int p = 0 ;
    int batch_size = 10000 ;
    lprint({"Extracting SFS strings on", to_string(config->threads), "threads.."});
    lprint({"Loading first batch.."});
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < config->threads; j++) {
            for (int k = 0; k <= batch_size / config->threads; k++) {
                bam_entries[i][j].push_back(bam_init1()) ;
            }
        }
    }
    load_batch_bam(config->threads, batch_size, p) ;
    batches.push_back(vector<vector<string>>(config->threads)) ; 
    //
    time_t t ;
    time(&t) ;
    int b = 0 ;
    uint64_t u = 0 ;

    bool should_load = true ;
    bool should_process = true ;
    bool loaded_last_batch = false ;
    bool should_update_current_batch = false ;

    uint64_t total_sfs = 0 ;
    uint64_t total_sfs_output_batch = 0 ;
    // load 
    current_input_batch = 0 ;
    load_input_sfs_batch() ;
    load_input_sfs_batch() ;
    assert(current_input_batch == sfs_batches.size()) ;

    while (true) {
        lprint({"Beginning batch", to_string(b + 1)});
        for (int i = 0 ; i < config->threads ; i++) {
            u += bam_entries[p][i].size() ;
        }
        if (!should_load) {
            should_process = false ;
        }
        if (loaded_last_batch) {
            should_load = false ;
        }
        #pragma omp parallel for num_threads(config->threads + 2)
        for(int i = 0; i < config->threads + 2; i++) {
            if (i == 0) {
                // load next batch of entries
                if (should_load) {
                    loaded_last_batch = !load_batch_bam(config->threads, batch_size, (p + 1) % 2) ;
                    lprint({"Loaded."});
                    batches.push_back(vector<vector<string>>(config->threads)) ; // previous and current output
                }
            } else if (i == 1) {
               if (b >= 1) {
                    // just count how many strings we have
                    uint64_t c = 0 ;
                    for (const auto &batch: batches[b - 1]) {
                        c += batch.size() ;
                    }
                    total_sfs += c ;
                    total_sfs_output_batch += c ;
                    cerr << "[I] Merged " << c << " new sequences. " << total_sfs << " total sequences. " << total_sfs_output_batch << " in current batch." << endl ;
                    // reached memory limit or last pipeline run
                    // output_batch(b - 1) ;
                }
            } else {
                // process current batch
                if (should_process) {
                    batches[b][i - 2] = process_batch(bam_entries[p][i - 2]) ; // mode==1 means input is bam: read sequence is already revcompled
                }
            }
        }

        if (total_sfs_output_batch >= 10000000 || !should_process) {
            output_batch(b) ;
            total_sfs_output_batch = 0 ;
            current_batch += 1 ;
            load_input_sfs_batch() ;
        }

        if (!should_load) {
            lprint({"Processed last batch of inputs."});
        }

        if (!should_process) {
            break ;
        }

        p += 1 ;
        p %= 2 ;
        b += 1 ;
        time_t s ;
        time(&s) ;
        if (s - t == 0) {
            s += 1 ;
        }
        cerr << "[I] Processed batch " << std::left << std::setw(10) << b << ". Reads so far " << std::right << std::setw(12) << u << ". Reads per second: " <<  u / (s - t) << ". Time: " << std::setw(8) << std::fixed << s - t << "\n" ;
    }
    sam_close(bam_file) ;
}

vector<string> Realigner::process_batch(vector<bam1_t*> &bam_entries) {
    vector<string> output ;
    int sfs_batch = -1 ;
    int buffer_len = 0 ;
    char* qseq = (char*) malloc(buffer_len) ;
    char* qqual = (char*) malloc(buffer_len) ;
    for (const auto alignment: bam_entries) { // a map
        string qname = string(bam_get_qname(alignment)) ;
        if (sfs_batch == -1) {
            if (sfs_batches[current_input_batch - 2].find(qname) != sfs_batches[current_input_batch - 2].end()) {
                sfs_batch = current_input_batch - 2 ;
            } else if (sfs_batches[current_input_batch - 1].find(qname) != sfs_batches[current_input_batch - 1].end()) {
                sfs_batch = current_input_batch - 1 ;
            } else {
                // this is a read that doesn't produce any SFS.
                continue ;
            }
        }
        // read data
        char *chrom = bam_header->target_name[alignment->core.tid];
        vector<pair<int, int>> alpairs = get_aligned_pairs(alignment) ;
        bool is_rev = bam_is_rev(alignment) ;
        uint flag = is_rev ? 16 : 0 ;
        uint32_t l = alignment->core.l_qseq ; //length of the read
        if (l >= buffer_len) {
            free(qseq) ;
            free(qqual) ;
            buffer_len = l + 1 ;
            qseq = (char*) malloc(buffer_len) ;
            qqual = (char*) malloc(buffer_len) ;
        }
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < l; i++){
            qseq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
            qqual[i] = ']'; // char(int(q[i]) + 33); //TODO: actual qual is lost
        }
        qseq[l] = '\0' ; // null terminate
        qqual[l] = '\0';
        for (const auto &sfs: sfs_batches[sfs_batch][qname]) {
            //cout << "processing SFS [" << sfs.s << ", " << sfs.s + sfs.l << "] on  read " << qname << endl ;
            vector<pair<int, int>> local_alpairs ;
            int start_pair_index = -1 ;
            int end_pair_index = 0 ;
            int j = 0 ;
            for (const auto &pos: alpairs) {
                if (pos.first != -1 && sfs.s <= pos.first && pos.first < sfs.s + sfs.l) {
                    local_alpairs.push_back(make_pair(pos.first, pos.second));
                    if (start_pair_index != -1) {
                        start_pair_index = j ;
                    }
                    end_pair_index = j ;
                }
                j++ ;
            }
            // TODO do we need this if?
            if (local_alpairs.empty()) {
                cerr << "EMPTY LOCAL " << qname << " " << sfs.s << " " << sfs.l << endl;
                continue;
            }
            // FILLING STARTING/ENDING -1s:
            // - if clips, we just add pairs til read end
            // - if insertion, we add pairs til first M we can find
            // slightly faster than looping over all pairs (can be thousands of extra iterations)
            if (local_alpairs.front().second == -1) {
                for (int j = start_pair_index - 1; j >= 0; j--) {
                    local_alpairs.insert(local_alpairs.begin(), alpairs[j]) ;
                    if (alpairs[j].second != -1) {
                        break ;
                    }
                }
            }
            //if (local_alpairs.front().second == -1) {
            //    bool add = false;
            //    for (int i = alpairs.size() - 1; i >= 0; --i) {
            //        if (add) { 
            //            local_alpairs.insert(local_alpairs.begin(), alpairs.at(i));
            //        }
            //        if (alpairs.at(i).second != -1 && add) {
            //            break;
            //        }
            //        if (!add && alpairs.at(i).first == local_alpairs.front().first) {
            //            add = true;
            //        }
            //    }
            //}
            if (local_alpairs.back().second == -1) {
                for (int j = end_pair_index + 1; j < alpairs.size(); j++) {
                    local_alpairs.push_back(alpairs[j]) ;
                    if (alpairs[j].second != -1) {
                        break ;
                    }
                }
            }
            //if (local_alpairs.back().second == -1) {
            //    bool add = false;
            //    for (int i = 0; i < alpairs.size(); ++i) {
            //        if (add) {
            //            local_alpairs.push_back(alpairs.at(i));
            //        }
            //        if (alpairs.at(i).second != -1 && add) {
            //            break;
            //        }
            //        if (!add && alpairs.at(i).first == local_alpairs.back().first) {
            //            add = true;
            //        }
            //    }
            //}
            uint qs = local_alpairs.front().first ;
            uint qe = local_alpairs.back().first ;
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
            if (local_alpairs.front().first == 0 && local_alpairs.front().second == -1) {
                // if we have initial clips, we get the position from original
                // alignments
                ts = alignment->core.pos;
            }
            if (local_alpairs.front().second == -1 && local_alpairs.back().second == -1) {
                cerr << "FULL CLIP " << qname << "." << sfs.s << ":" << sfs.l << endl;
                continue;
            }
            //
            cout << chrom << endl ;
            CIGAR localcigar = rebuild_cigar(chromosome_seqs[chrom], qseq, local_alpairs);
            localcigar.fixclips();
            string localqseq(qseq + qs, qe - qs + 1);
            string localqqual(qqual + qs, qe - qs + 1);
            // maybe do something more elegant
            string o = qname + "." + to_string(sfs.s) + "-" + to_string(sfs.s + sfs.l - 1) + "\t"
                + to_string(flag) + "\t"
                + chrom + "\t"
                + to_string(ts + 1) + "\t"
                + to_string(alignment->core.qual) + "\t"
                + localcigar.to_str() + "\t"
                + "*" + "\t"
                + "0" + "\t"
                + "0" + "\t"
                + localqseq + "\t"
                + localqqual + "\t"
                + "NM:i:" + to_string(localcigar.mismatches) ;
            output.push_back(o) ;
        }
    }
    free(qseq) ;
    free(qqual) ;
    return output ; 
}

void Realigner::output_batch(int b) {
    for (int j = last_dumped_batch; j < b; j++) {
        for (int i = 0; i < config->threads; i++) {
            for (auto f: batches[j][i]) {
                out_file << f << endl ;
            }
        }
        batches[j].clear() ;
        last_dumped_batch++ ;
    }
}

CIGAR Realigner::rebuild_cigar(const string &ref_seq, const string &read_seq, const vector<pair<int, int>> &alpairs) {
    CIGAR cigar;
    int last_ref_pos = alpairs[0].second;
    for (const pair<int, int> p: alpairs) {
        int ref_pos = p.second ;
        int read_pos = p.first ;
        if (last_ref_pos != -1 && ref_pos != -1) {
            if (last_ref_pos != ref_pos && last_ref_pos != ref_pos - 1) {
                // Deletions
                int d = ref_pos - last_ref_pos - 1;
                cigar.add(d, 'D', d);
            }
        }
        if (read_pos != -1 && ref_pos != -1) {
            // (mis)Match
            if (ref_seq[ref_pos] != read_seq[read_pos]) {
                cigar.add(1, 'M', 1);
            } else {
                cigar.add(1, 'M', 0);
            }
        } else if (ref_pos == -1) {
            // Insertion (1bp)
            cigar.add(1, 'I', 1) ;
        }
        last_ref_pos = ref_pos ;
    }
    return cigar;
}
