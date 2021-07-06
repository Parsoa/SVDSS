#ifndef TAU_HPP
#define TAU_HPP

#include <omp.h>
#include <vector>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <pthread.h>
#include <unordered_map>

#include "rle.h"
#include "rld0.h"
#include "mrope.h"

#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"

#include "fastq.hpp"
#include "lprint.hpp"
#include "config.hpp"
#include "bed_utils.hpp"

using namespace std;

struct TauDistribution {
    std::unordered_map<string, std::vector<Locus>> loci ;
    std::vector<std::string> sfsnames ;
    int num_loci ;
};

class Tau {

public:
    void run();

private:

    int current_batch = 0 ;

    std::unordered_map<std::string, TauDistribution> tau_distribution ;

    Configuration* config ;
};

#endif
