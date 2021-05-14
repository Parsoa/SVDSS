#ifndef CNT_HPP
#define CNT_HPP

#include <string>
#include <unordered_map>
#include "omp.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "config.hpp"
#include "lprint.hpp"

using namespace std ;

class Counter {

public:

    void load_counts(std::string path) ;

    std::unordered_map<std::string, int> counts ;

} ;

#endif
