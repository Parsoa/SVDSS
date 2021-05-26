#ifndef ASSEMBLER_HPP
#define ASSEMBLER_HPP

#include "config.hpp"
#include "sfsutils.hpp"
#include <fstream>
#include <list>

using namespace std;

class Assembler {
public:
  void run();
  vector<SFS> assemble(vector<SFS> &sfs);
};

#endif
