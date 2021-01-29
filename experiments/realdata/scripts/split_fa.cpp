// g++ -std=c++11 -O3 split.cpp -o split -lz

#include <iostream>
#include <fstream>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;

int main(int argc, char *argv[])
{
    char *fq_path = argv[1];
    int l1 = 500; // atoi(argv[2]);
    // int l1 = 70; // atoi(argv[2]);
    // int l2 = 1000; // atoi(argv[3]);
    string prefix_path = argv[2];

    ofstream shortspec(prefix_path + ".short.fa");
    // ofstream midspec (prefix_path + ".mid.fq");
    ofstream longspec(prefix_path + ".long.fa");

    gzFile fp = gzopen(fq_path, "r");
    kseq_t *seq = kseq_init(fp);
    int l;
    while ((l = kseq_read(seq)) >= 0)
    {
        if (l <= l1)
            shortspec << ">" << seq->name.s << "\n"
                      << seq->seq.s << "\n";
        // else if (l1 <= l && l <= l2)
        //   midspec << "@" << seq->name.s << "\n"
        //           << seq->seq.s << "\n"
        //           << "+"
        //           << "\n"
        //           << seq->qual.s << "\n";
        else
            longspec << ">" << seq->name.s << "\n"
                     << seq->seq.s << "\n";
    }
    kseq_destroy(seq);
    gzclose(fp);
    shortspec.close();
    // midspec.close();
    longspec.close();
    return 0;
}
