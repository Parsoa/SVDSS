#include "sfs.hpp"

using namespace std;

bool operator<(const SFS &x, const SFS &y) { return x.s < y.s; }

map<string, vector<SFS>> parse_sfsfile(const string &sfs_path, int tau) {
  map<string, vector<SFS>> SFSs;
  string line;
  ifstream inf(sfs_path);
  if (inf.is_open()) {
    string info[5];
    string read_name;
    while (getline(inf, line)) {
      stringstream ssin(line);
      int i = 0;
      while (ssin.good() && i < 5) {
        ssin >> info[i++];
      }
      if (stoi(info[3]) < tau) {
        continue;
      }
      if (info[0].compare("*") != 0) {
        read_name = info[0];
        SFSs[read_name] = vector<SFS>();
      }
      // TODO
      SFSs[read_name].push_back(SFS(stoi(info[1]), stoi(info[2]), stoi(info[3]),
                                    true)); // info[4].compare("1") == 0));
    }
  }
  return SFSs;
}
