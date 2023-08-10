#include "sfs.hpp"

// TODO: if we split the .sfs into more files, we can load it using more threads
// (but for SVs having a single file shouldn't be the bottleneck)
unordered_map<string, vector<SFS>> parse_sfsfile(const string &sfs_path) {
  spdlog::info("Loading SFSs from {}..", sfs_path);
  unordered_map<string, vector<SFS>> SFSs;
  int total = 0;
  string line;
  ifstream inf(sfs_path);
  if (inf.is_open()) {
    string info[4];
    string read_name;
    while (getline(inf, line)) {
      stringstream ssin(line);
      int i = 0;
      while (ssin.good() && i < 4)
        ssin >> info[i++];
      if (info[0].compare("*") != 0) {
        read_name = info[0];
        SFSs[read_name] = vector<SFS>();
      }
      SFSs[read_name].push_back(
          SFS(read_name, stoi(info[1]), stoi(info[2]), stoi(info[3])));
      ++total;
    }
  }
  spdlog::info("Loaded {} SFSs from {} reads.", total, SFSs.size());
  return SFSs;
}
