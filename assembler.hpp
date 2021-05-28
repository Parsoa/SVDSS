#ifndef ASSEMBLER_HPP
#define ASSEMBLER_HPP

#include <list>
#include <fstream>

#include "config.hpp"
#include "sfsutils.hpp"

using namespace std;

class Assembler {
public:
  void run();
  vector<SFS> assemble(vector<SFS> &sfs);
};

#endif
