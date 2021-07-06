#include "realigner.hpp"
#include "chromosomes.hpp"

bool Realigner::load_batch_bam(int threads, int batch_size, int p) {
    int i = 0 ;
    int n = 0 ;
    int m = 0 ;
    while (sam_read1(bam_file, bam_header, bam_entries[p][n % threads][i]) >= 0) {
        m++ ;
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
        string qname = string(bam_get_qname(alignment)) ;
        if (SFSs.find(qname) == SFSs.end()) {
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
    lprint({"Loaded", to_string(n), "BAM reads with SFS from total of", to_string(m), "reads."});
    return n != 0 ? true : false ;
}

void Realigner::load_input_sfs_batch() {
    int threads = config->threads ; 
    int num_batches = config->aggregate_batches ;
    int num_threads = num_batches < threads ? num_batches : threads ;
    vector<unordered_map<string, vector<SFS>>> _SFSs(num_batches) ;
    cout << "Loading assmbled SFS.." << endl ;
    #pragma omp parallel for num_threads(num_threads)
    for (int j = 0; j < num_batches; j++) {
        string s_j = std::to_string(j) ;
        string inpath = config->workdir + "/solution_batch_" + s_j + ".assembled.sfs" ;
        cout << "[I] Loading SFS from " << inpath << endl ;
        ifstream inf(inpath) ;
        string line ;
        if (inf.is_open()) {
            string info[4];
            string read_name;
            while (getline(inf, line)) {
                stringstream ssin(line);
                int i = 0;
                while (ssin.good() && i < 4) {
                    ssin >> info[i++];
                }
                if (info[0].compare("*") != 0) {
                    read_name = info[0];
                    _SFSs[j][read_name] = vector<SFS>();
                }
                _SFSs[j][read_name].push_back(SFS(stoi(info[1]), stoi(info[2]), stoi(info[3]), true)) ;
            }
        }
    }
    int r = 0 ;
    int c = 0 ;
    for (int j = 0; j < num_batches; j++) {
        lprint({"Batch", to_string(j), "with", to_string(_SFSs[j].size()), "strings."});
        r += _SFSs[j].size() ;
        SFSs.insert(_SFSs[j].begin(), _SFSs[j].end()) ;
        for (auto& read: _SFSs[j]) {
            c += read.second.size() ;
        }
    }
    lprint({"Aligning", to_string(c), "SFS strings on", to_string(r), " reads.", to_string(config->threads), "threads.."}) ;
}

void Realigner::load_target_SFS_set() {
    ifstream txt_file("num_loci.txt") ;
    string line ;
    while (getline(txt_file, line)) {
        int p = line.rfind(':') ;
        int count = std::stoi(line.substr(p + 1, line.length() - (p + 1))) ;
        string seq = line.substr(0, p) ;
        if (count > 1) {
            target_sfs[seq] = count ;
        }
    }
    cout << "Loaded " << target_sfs.size() << " target SFS sequences." << endl ;
}

void Realigner::run() {
    config = Configuration::getInstance() ;
    if (config->target != "") {
        load_target_SFS_set() ;
    }
    load_input_sfs_batch() ; // load all SFSs
    load_chromosomes(config->reference) ;
    out_file = ofstream(config->workdir + "/realignments.sam") ;
    tau_out_file = ofstream(config->workdir + "/superstring_loci.bed") ;
    bam_file = sam_open(config->bam.c_str(), "r") ;
    bam_header = sam_hdr_read(bam_file) ; //read header
    out_file << bam_header->text; // sam_hdr_str(bamhdr) "was not declared in this scope"
    int batch_size = 10000 ;
    for (int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(config->threads)) ;
        for (int j = 0; j < config->threads; j++) {
            for (int k = 0; k < batch_size / config->threads; k++) {
                bam_entries[i][j].push_back(bam_init1()) ;
            }
        }
    }
    lprint({"Loading first batch.."});
    batches.push_back(vector<vector<batch_atom_type>>(config->threads)) ; 
    //
    time_t t ;
    time(&t) ;
    int p = 0 ;
    int b = 0 ;
    uint64_t u = 0 ;

    bool should_load = true ;
    bool should_process = true ;
    bool loaded_last_batch = false ;
    bool should_update_current_batch = false ;

    uint64_t total_sfs = 0 ;
    uint64_t total_sfs_output_batch = 0 ;
    // load 
    load_batch_bam(config->threads, batch_size, p) ;
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
                    batches.push_back(vector<vector<batch_atom_type>>(config->threads)) ; // previous and current output
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
                    cerr << "[I] Merged " << c << " new alignments. " << total_sfs << " total alignments. " << total_sfs_output_batch << " in current batch." << endl ;
                }
            } else {
                // process current batch
                if (should_process) {
                    batches[b][i - 2] = process_batch(p, i - 2);
                }
            }
        }

        if (total_sfs_output_batch >= 10000000 || !should_process) {
            lprint({"Memory limit reached, dumping alignments.."}) ;
            output_batch(b) ;
            total_sfs_output_batch = 0 ;
            current_batch += 1 ;
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

vector<batch_atom_type> Realigner::process_batch(int p, int index) {
    vector<batch_atom_type> output ;
    int buffer_len = 20000 ;
    char* qseq = (char*) malloc(buffer_len) ;
    char* qqual = (char*) malloc(buffer_len) ;
    for (int b = 0; b < bam_entries[p][index].size(); b++) {
        bam1_t* alignment = bam_entries[p][index][b] ;
        string qname = string(bam_get_qname(alignment)) ;
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
            qqual[i] = ']'; // char(int(q[i]) + 33);
        }
        qseq[l] = '\0' ; // null terminate
        qqual[l] = '\0';
        int last_pos = 0 ;
        for (int s = 0; s < SFSs[qname].size(); s++) {
            auto& sfs = SFSs[qname][s] ;
            //
            vector<pair<int, int>> local_alpairs ;
            int start_pair_index = -1 ;
            int end_pair_index = -1 ;
            for (int i = last_pos; i < alpairs.size(); i++) {
                if (alpairs[i].first != -1 && alpairs[i].first >= sfs.s && alpairs[i].first < sfs.s + sfs.l) {
                    local_alpairs.push_back(make_pair(alpairs[i].first, alpairs[i].second));
                    if (start_pair_index == -1) {
                        start_pair_index = i ;
                    }
                    end_pair_index = i ;
                }
                if (alpairs[i].first != -1 && alpairs[i].first > sfs.s + sfs.l) {
                    break ;
                }
            }
            //
            string sfsseq(qseq + sfs.s, sfs.l) ;
            string sfsqual(qqual + sfs.s, sfs.l) ;
            if (config->target != "") {
                if (target_sfs.find(sfsseq) == target_sfs.end()) {
                    continue ;
                }
            }

            if (local_alpairs.empty()) {
                assert(end_pair_index == -1 && start_pair_index == -1) ;
                //cerr << "EMPTY LOCAL " << qname << " " << sfs.s << " " << sfs.l << endl;
                continue ;
            }
            if (s < SFSs[qname].size() - 1) {
                auto& next = SFSs[qname][s + 1] ;
                if (next.s >= sfs.s + sfs.l) {
                    last_pos = end_pair_index - 1 ;
                } else {
                    last_pos = start_pair_index ;
                }
            }
            // FILLING STARTING/ENDING -1s:
            // - if clips, we just add pairs til read end
            // - if insertion, we add pairs til first M we can find
            if (local_alpairs[0].second == -1) {
                for (int j = start_pair_index - 1; j >= 0; j--) {
                    local_alpairs.insert(local_alpairs.begin(), alpairs[j]) ;
                    if (alpairs[j].second != -1) {
                        break ;
                    }
                }
            }
            if (local_alpairs.back().second == -1) {
                for (int j = end_pair_index + 1; j < alpairs.size(); j++) {
                    local_alpairs.push_back(alpairs[j]) ;
                    if (alpairs[j].second != -1) {
                        break ;
                    }
                }
            }

            uint qs = local_alpairs.front().first ;
            uint qe = local_alpairs.back().first ;
            // In some (very rare I hope) cases, an insertion follows a deletions (or
            // viceversa). So we are trying to find the first M - that is a non -1 in
            // the pairs - but that pair has -1 on the read
            if (qs == -1 || qe == -1) {
                //cerr << "INS-DEL " << qname << "." << sfs.s << ":" << sfs.l << endl;
                continue;
            }
            // If clips, we have trailing -1 in target positions. We have to find the first placed base
            int ts = local_alpairs.front().second;
            if (local_alpairs.front().first == 0 && local_alpairs.front().second == -1) {
                // if we have initial clips, we get the position from original alignments
                ts = alignment->core.pos;
            }
            if (local_alpairs.front().second == -1 && local_alpairs.back().second == -1) {
                //cerr << "FULL CLIP " << qname << "." << sfs.s << ":" << sfs.l << endl;
                continue;
            }
            CIGAR localcigar = rebuild_cigar(chromosome_seqs[chrom], qseq, local_alpairs);
            localcigar.fixclips();
            string localqseq(qseq + qs, qe - qs + 1);
            string localqqual(qqual + qs, qe - qs + 1);
            string sfsname = qname + "." + to_string(sfs.s) + "-" + to_string(sfs.s + sfs.l - 1) ;
            string o = sfsname + "\t"
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
            string tau_o = sfsname + "\t" + sfsseq + "\t" + chrom + "\t" + to_string(ts + 1) ;
            output.push_back(make_pair(o, tau_o)) ;
        }
    }
    free(qseq) ;
    free(qqual) ;
    return output ; 
}

void Realigner::output_batch(int b) {
    cout << "Outputting batch " << b << ".." << endl ;
    for (int j = last_dumped_batch; j < b; j++) {
        for (int i = 0; i < config->threads; i++) {
            for (auto f: batches[j][i]) {
                out_file << f.first << endl ;
                tau_out_file << f.second << endl ;
            }
        }
        batches[j].clear() ;
        last_dumped_batch++ ;
    }
}

CIGAR Realigner::rebuild_cigar(char* ref_seq, char* read_seq, const vector<pair<int, int>> &alpairs) {
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
