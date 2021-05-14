#include "counter.hpp"

void Counter::load_counts(string base) {
    auto c = Configuration::getInstance() ;
    cout << "loading counts.." << endl ;
    vector<unordered_map<string, int>> _counts(c->batch_size) ;
    #pragma omp parallel for num_threads(c->threads)
    for (int j = 0; j < c->batch_size; j++) {
        string s_j = std::to_string(j) ;
        string index = j < 10 ? "00" + s_j : j < 100 ? "0" + s_j : s_j ;
        string path = base + "/solution/sequence_counts_batch_" + index + ".txt" ;
        ifstream txt_file(path) ;
        string line ;
        while (std::getline(txt_file, line)) {
            stringstream ss(line) ;
            vector<string> tokens ;
            string token ;
            while (getline(ss, token, '\t')) {
                tokens.push_back(token) ;
            }
            int p = tokens[0].rfind(':') ;
            string name = tokens[0].substr(0, p - 1) ;
            if (_counts[j].find(name) == _counts[j].end()) {
                _counts[j][name] = 0 ;
            }
            _counts[j][name] += std::stoi(tokens[1]) ;
        }
        txt_file.close() ;
    }
    cout << "Merging.." << endl ;
    for (int j = 0; j < c->batch_size; j++) {
        for (auto it = _counts[j].begin(); it != _counts[j].end(); it++) {
            if (counts.find(it->first) == counts.end()) {
                counts[it->first] = 0 ;
            }
            counts[it->first] += it->second ;
        }
    }
    cout << "Loaded " << counts.size() << " counts." << endl ;
}
