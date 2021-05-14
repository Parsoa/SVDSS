#include "sfsutils.hpp"

bool operator<(const SFS &x, const SFS &y) { return x.s < y.s; }

map<string, vector<SFS>> parse_sfsfile(const string &sfs_path, int tau) {
    map<string, vector<SFS>> SFSs;
    string line;
    ifstream inf(sfs_path);
    if (inf.is_open()) {
        string info[4];
        while (getline(inf, line)) {
            stringstream ssin(line);
            int i = 0;
            while (ssin.good() && i < 4) {
                ssin >> info[i++];
            }
    
            if (stoi(info[3]) < tau) {
                continue;
            }
    
            if (SFSs.find(info[0]) == SFSs.end()) {
                SFSs[info[0]] = vector<SFS>();
            }
            SFSs[info[0]].push_back(SFS(stoi(info[1]), stoi(info[2]), stoi(info[3])));
        }
    }
    return SFSs ;
}
