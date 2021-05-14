#include "aggregator.hpp"

void Aggregator::run() {
    auto c = Configuration::getInstance() ;
    find_high_abundance_sequences() ;
    load_sequences() ;
    dump_sequences() ;
}

void Aggregator::find_high_abundance_sequences() {
    auto c = Configuration::getInstance() ;
    num_batches = c->aggregate_batches ;
    lprint({"Finding high-abundance strings from", to_string(num_batches), "batches.."});
    vector<unordered_map<int, int>> _sequences(num_batches) ;
    int e = 0 ;
    // cout first pass
    #pragma omp parallel for num_threads(num_batches)
    for (int j = 0; j < num_batches; j++) {
        string s_j = std::to_string(j) ;
        string path = c->workdir + "/solution_batch_" + s_j + ".fastq" ;
        ifstream txt_file(path) ;
        int i = 0 ;
        int count ;
        string line ;
        string name ;
        while (std::getline(txt_file, line)) {
            if (i == 0) {
                int p = line.rfind(':') ;
                count = std::stoi(line.substr(p + 1, line.length() - (p + 1))) ;
            }
            if (i == 1) {
                string canon = canonicalize(line) ;
                int hash = std::hash<std::string>()(canon) ;
                if (_sequences[j].find(hash) == _sequences[j].end()) {
                    _sequences[j][hash] = 0 ;
                }
                _sequences[j][hash] += count ;
            }
            i++ ;
            i %= 4 ;
        }
        txt_file.close() ;
    }
    lprint({"Merging batches.."});
    for (int j = 0; j < num_batches; j++) {
      lprint({"Batch", to_string(j), "with", to_string(_sequences[j].size()), "sequences."});
        for (auto it = _sequences[j].begin(); it != _sequences[j].end(); it++) {
            if (sequence_index.find(it->first) == sequence_index.end()) {
                sequence_index[it->first] = 0 ;
            }
            sequence_index[it->first] += it->second ;
        }
        _sequences[j].clear() ;
    }
    lprint({"Loaded", to_string(sequence_index.size()), "SFS sequences."});
    auto it = sequence_index.begin();
    while (it != sequence_index.end()) {
        if (it->second < c->cutoff) {
            it = sequence_index.erase(it) ; 
        } else {
            it++ ;
        }
    }
    lprint({to_string(sequence_index.size()), "high-abundance sequences found."});;
}

void Aggregator::load_sequences() {
    auto c = Configuration::getInstance() ;
    num_batches = c->aggregate_batches ;
    lprint({"Loading sequences from", to_string(num_batches), "batches.."});
    vector<unordered_map<string, int>> _sequences(num_batches) ;
    vector<unordered_map<string, unordered_map<string, int>>> _read_ids(num_batches) ;
    int e = 0 ;
    // cout first pass
    #pragma omp parallel for num_threads(num_batches)
    for (int j = 0; j < num_batches; j++) {
        string s_j = std::to_string(j) ;
        string path = c->workdir + "/solution_batch_" + s_j + ".fastq" ;
        ifstream txt_file(path) ;
        int i = 0 ;
        int count ;
        string line ;
        string name ;
        while (std::getline(txt_file, line)) {
            if (i == 0) {
                int p = line.rfind(':') ;
                name = line.substr(0, p - 1) ;
                count = std::stoi(line.substr(p + 1, line.length() - (p + 1))) ;
            }
            if (i == 1) {
                string canon = canonicalize(line) ;
                int hash = std::hash<std::string>()(canon) ;
                if (sequence_index.find(hash) != sequence_index.end()) {
                    if (_sequences[j].find(canon) == _sequences[j].end()) {
                        _sequences[j][canon] = 0 ;
                    }
                    _sequences[j][canon] += count ;
                }
            }
            i++ ;
            i %= 4 ;
        }
        // load read ids
        string id_path = c->workdir + "/read_ids_batch_" + s_j + ".fasta" ;
        ifstream id_txt_file(id_path) ;
        i = 0 ;
        string canon ;
        while (std::getline(id_txt_file, line)) {
            if (i == 0) {
                canon = canonicalize(line.substr(1, line.size() - 1)) ;
                //assert(_sequences[j].find(canon) != _sequences[j].end()) ;
            }
            if (i == 1) {
                int hash = std::hash<std::string>()(canon) ;
                if (sequence_index.find(hash) != sequence_index.end()) {
                    _read_ids[j][canon][line] += 1 ;
                }
            }
            i++ ;
            i %= 2 ;
        }
        txt_file.close() ;
        id_txt_file.close() ;
    }
    lprint({"Merging.."});
    for (int j = 0; j < num_batches; j++) {
      lprint({"Batch", to_string(j), "with", to_string(_sequences[j].size()), "sequences."});
        for (auto it = _sequences[j].begin(); it != _sequences[j].end(); it++) {
            if (sequences.find(it->first) == sequences.end()) {
                sequences[it->first] = 0 ;
            }
            sequences[it->first] += it->second ;
            read_ids[it->first].insert(_read_ids[j][it->first].begin(), _read_ids[j][it->first].end()) ;
        }
        _sequences[j].clear() ;
    }
    lprint({"Loaded", to_string(sequences.size()), "seuqences."});
}

void Aggregator::dump_sequences() {
    auto c = Configuration::getInstance() ;
    lprint({"Dumping aggregated counts.."});
    string path = c->workdir + "/solution_aggregated.fastq" ;
    string id_path = c->workdir + "/read_ids_aggregated.fasta" ;
    ofstream o(path) ;
    ofstream id_o(id_path) ;
    int i = 0 ;
    int n = 0 ;
    for (auto it = sequences.begin(); it != sequences.end(); it++) {
        if (it->second >= c->cutoff) {
            o << "@sol_" << i << "#" << it->second << endl ;
            o << it->first << endl ;
            o << "+" << endl ;
            for (int j = 0; j < it->first.length(); j++) {
                o << "I" ;
            }
            o << endl ;
            // read ids
            id_o << it->first << endl ;
            for (auto s: read_ids[it->first]) {
                id_o << s.first ;
            }
            id_o << endl ;
            //
            n += 1 ;
        }
        i += 1 ;
    }
    lprint({to_string(n), "useful sequences remaining."});
}
