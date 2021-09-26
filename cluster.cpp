#include "cluster.hpp"

using namespace std ;

unsigned char _char26_table[256] = {
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 0, 5, 1, 6, 7, 8, 2, 9, 10, 11, 12, 13, 14, 4, 15,
    16, 17, 18, 19, 3, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26,
    26, 0, 5, 1, 6, 7, 8, 2, 9, 10, 11, 12, 13, 14, 4, 15,
    16, 17, 18, 19, 3, 20, 21, 22, 23, 24, 25, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
    26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26};

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
    abpt->is_diploid = 1;
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
    uint8_t **cons_seq;
    int **cons_cov, *cons_l, cons_n = 0;
    uint8_t **msa_seq;
    int msa_l = 0;
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);

    string cons = "";
    for (int j = 0; j < cons_l[0]; ++j) {
        cons += "ACGTN-"[cons_seq[0][j]] ;
    }
    if (cons_n > 0) {
        for (int i = 0; i < cons_n; ++i) {
            free(cons_seq[i]);
        }
        free(cons_seq);
        // free(cons_cov);
        free(cons_l);
    }
    if (msa_l) {
        for (uint i = 0; i < n_seqs; ++i) {
            free(msa_seq[i]);
        }
        free(msa_seq);
    }
    for (uint i = 0; i < n_seqs; ++i) {
        free(bseqs[i]);
    }
    free(bseqs);
    free(seq_lens);
    abpoa_free(ab);
    abpoa_free_para(abpt);
    return cons ;
}

void Cluster::dump(ofstream& o) const {
    //o << chrom << " " << s << " " << e << " " << cov << " " << size();
    for (const string &seq: seqs) {
        o << " " << seq;
    }
    o << endl ;
}
