#include "converter.hpp"

bool Converter::load_batch_bam(int threads, int batch_size, int p) {
    for (int i = 0; i < threads; i++) {
        fastq_entries[p][i].clear() ;
    }
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
        uint32_t l = alignment->core.l_qseq ; //length of the read
        char* seq = (char*) malloc(l + 1) ;
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < l; i++){
            seq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        seq[l] = '\0' ; // null terminate
        fastq_entries[p][n % threads].push_back(fastq_entry_t(bam_get_qname(alignment), seq, seq)) ;
        n += 1 ;
        if (n % threads == threads - 1) {
            i += 1 ;
        }
        if (n == batch_size) {
            return true ;
        }
    }
    lprint({"Loaded", to_string(n), "BAM reads.."});
    return n != 0 ? true : false ;
}

bool Converter::load_batch_fastq(int threads, int batch_size, int p) {
    for (int i = 0; i < threads; i++) {
        fastq_entries[p][i].clear() ;
    }
    int l = 0 ;
    int n = 0 ;
    while ((l = kseq_read(fastq_iterator)) >= 0) {
        fastq_entries[p][n % threads].push_back(fastq_entry_t(fastq_iterator->name.s, fastq_iterator->seq.s, fastq_iterator->qual.s)) ;
        n += 1 ;
        if (n == batch_size) {
          lprint({"Loaded", to_string(n), "FASTQ reads.."});
            return true ;
        }
    }
    lprint({"Loaded", to_string(n), "FASTQ reads.."});
    // TODO: this results in an extra empty batch if total number of reads is not a multiple of batch size
    return n != 0 ? true : false ;
}

// load SFS file and round-robin reads between threads
void Converter::load_input_sfs_batch() {
    if (current_input_batch == config->aggregate_batches) {
        cerr << "[I] all input SFS batches loaded. Not loading new batch." << endl ;
        return ;
    }
    string s_j = std::to_string(current_input_batch) ;
    string inpath = config->workdir + "/solution_batch_" + s_j + ".sfs" ;
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

void Converter::run() {
    // Load all SFS in current input batch
    // Load reads and round-robin between threads
    config = Configuration::getInstance() ;
    // parse arguments
    lprint({"Done."});
    int mode = 0 ;
    if (config->fastq != "") {
        lprint({"FASTQ input:", config->fastq}) ;
        fastq_file = gzopen(config->fastq.c_str(), "r") ;
        fastq_iterator = kseq_init(fastq_file) ;
    } else if (config->bam != "") {
        lprint({"BAM input.."});
        bam_file = hts_open(config->bam.c_str(), "r") ;
        bam_header = sam_hdr_read(bam_file) ; //read header
        mode = 1 ;
    } else {
        lprint({"No input file provided, aborting.."}, 2);
        exit(1) ;
    }
    // load first batch
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(config->threads)) ;
        fastq_entries.push_back(vector<vector<fastq_entry_t>>(config->threads)) ; // current and next output
    }
    int p = 0 ;
    int batch_size = 10000 ;
    lprint({"Extracting SFS strings on", to_string(config->threads), "threads.."});
    lprint({"Loading first batch.."});
    if (mode == 0) {
        load_batch_fastq(config->threads, batch_size, p) ;
    } else {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < config->threads; j++) {
                for (int k = 0; k <= batch_size / config->threads; k++) {
                    bam_entries[i][j].push_back(bam_init1()) ;
                }
            }
        }
        load_batch_bam(config->threads, batch_size, p) ;
    }
    batches.push_back(vector<vector<fastq_entry_t>>(config->threads)) ; 
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
            u += fastq_entries[p][i].size() ;
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
                    if (mode == 0) {
                        loaded_last_batch = !load_batch_fastq(config->threads, batch_size, (p + 1) % 2) ;
                    } else {
                        loaded_last_batch = !load_batch_bam(config->threads, batch_size, (p + 1) % 2) ;
                    }
                    lprint({"Loaded."});
                    batches.push_back(vector<vector<fastq_entry_t>>(config->threads)) ; // previous and current output
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
                }
            } else {
                // process current batch
                if (should_process) {
                    batches[b][i - 2] = process_batch(fastq_entries[p][i - 2]) ; // mode==1 means input is bam: read sequence is already revcompled
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

    // cleanup
    kseq_destroy(fastq_iterator) ;
    gzclose(fastq_file) ;
}

vector<fastq_entry_t> Converter::process_batch(vector<fastq_entry_t> &fastq_entries) {
    vector<fastq_entry_t> output ;
    int sfs_batch = -1 ;
    for (const auto &read: fastq_entries) { // a map
        if (sfs_batch == -1) {
            if (sfs_batches[current_input_batch - 2].find(read.head) != sfs_batches[current_input_batch - 2].end()) {
                sfs_batch = current_input_batch - 2 ;
            } else if (sfs_batches[current_input_batch - 1].find(read.head) != sfs_batches[current_input_batch - 1].end()) {
                sfs_batch = current_input_batch - 1 ;
            } else {
                // this is a read that doesn't produce any SFS.
                continue ;
            }
        }
        for (const auto &sfs: sfs_batches[sfs_batch][read.head]) {
            string header = read.head + "#" + to_string(sfs.s) + "#" + to_string(sfs.s + sfs.l - 1) + "#" + to_string(sfs.c) ;
            string sfsseq(read.seq, sfs.s, sfs.l) ;
            string sfsqual(read.qual, sfs.s, sfs.l) ;
            output.push_back(fastq_entry_t(header, sfsseq, sfsqual)) ; 
        }
    }
   return output ; 
}

void Converter::output_batch(int b) {
    auto c = Configuration::getInstance();
    string path = c->workdir + "/solution_batch_" + std::to_string(current_batch) + ".fastq";
    lprint({"Outputting to", path}) ;
    std::ofstream o(path) ;
    for (int j = last_dumped_batch; j < b; j++) {
        for (int i = 0; i < config->threads; i++) {
            for (auto f: batches[j][i]) {
                o << "#" + f.head << endl ;
                o << f.seq << endl ;
                o << "+" << endl ;
                o << f.qual << endl ;
            }
        }
        batches[j].clear() ;
        last_dumped_batch++ ;
    }
}
