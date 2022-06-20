#include "cluster.hpp"

using namespace std ;

// AaCcGgTtNn ==> 0,1,2,3,4
unsigned char _char26_table[256] = {
       0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
       4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

string Cluster::poa() {
    // --- POA
    uint n_seqs = seqs.size();
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->disable_seeding = 1;
    abpt->align_mode = 0; // global
    abpt->out_msa = 0;
    abpt->out_cons = 1;
    abpt->out_gfa = 0;
    // abpt->w = 6, abpt->k = 9;
    // abpt->min_w = 10; // minimizer-based seeding and partition
    // abpt->is_diploid = 1;
    // abpt->min_freq = 0.3;
    abpt->progressive_poa = 0;
    abpoa_post_set_para(abpt);

    int *seq_lens = (int *)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * n_seqs);
    for (uint i = 0; i < n_seqs; ++i) {
        seq_lens[i] = seqs[i].size();
        bseqs[i] = (uint8_t *)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (int j = 0; j < seq_lens[i]; ++j) {
            bseqs[i][j] = _char26_table[(int)seqs[i][j]];
        }
    }

    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, NULL);

    abpoa_cons_t *abc = ab->abc;

    string cons = "";
    if (abc->n_cons > 0) {
        for (int j = 0; j < abc->cons_len[0]; ++j) {
            cons += "ACGTN"[abc->cons_base[0][j]];
        }
    }

    for (uint i = 0; i < n_seqs; ++i) { free(bseqs[i]); }
    free(bseqs); free(seq_lens);
    abpoa_free(ab); abpoa_free_para(abpt);

    return cons ;
}

void Cluster::dump(ofstream& o) const {
    //o << chrom << " " << s << " " << e << " " << cov << " " << size();
    for (const string &seq: seqs) {
        o << " " << seq;
    }
    o << endl ;
}
