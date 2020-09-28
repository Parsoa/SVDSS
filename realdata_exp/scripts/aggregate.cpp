#include <iostream>
#include <fstream>

#include <unordered_map>

using namespace std;

int main(int argc, char *argv[]) {
  string dir_path = argv[1];
  unordered_map<string, int> seqs;
  
  ifstream in_file;
  string line;
  int i, nl, p, count;
  for(i=0; i<2048; ++i) { // hardcoded
    in_file.open(dir_path + "/solution_batch_" + to_string(i) + ".fastq");
    if(!in_file.is_open()) {
      in_file.close();
      break;
    }
    nl = 0;
    while(getline(in_file, line)) {
      if (nl == 0) {
	p = line.rfind(':');
	count = stoi(line.substr(p + 1, line.length() - (p + 1)));
      } else if(nl == 1) {
	if (seqs.find(line) == seqs.end())
	  seqs[line] == 0;
	seqs[line] += count ;
      }
      ++nl;
      nl %= 4;
    }
    in_file.close();
  }

  i = 0;
  for(auto it = seqs.begin(); it != seqs.end(); ++it) {
    if (it->second > 1) {
      cout << "@sol_" << i << '#' << it->second << endl ;
      cout << it->first << endl ;
      cout << "+" << endl ;
      string qual (it->first.size(), ']');
      cout << qual << endl;
      ++i;
    }
  }
}
