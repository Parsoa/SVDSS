#ifndef ASSEMBLER_HPP
#define ASSEMBLER_HPP

#include <fstream>
#include <list>

#include "config.hpp"
#include "sfs.hpp"

using namespace std;

class Assembler {
public:
  void run();
  vector<SFS> assemble(vector<SFS> &sfs);
};

#endif
