#include "cluster.hpp"

using namespace std ;

const unsigned char _nt4_table[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

Cluster::Cluster(const string &chrom_) {
    chrom = chrom_;
    s = -1;
    e = 0;
}

void Cluster::add_fragment(Fragment f) {
    // TODO improve this
    //if (find(fragments.begin(), fragments.end(), f) == fragments.end()) {
    fragments.push_back(f) ;
    s = min(s, f.ref_s) ;
    e = max(e, f.ref_e) ;
    //}
}

void Cluster::set_full_coverage(const uint c) {
    full_cov = c;
}

string Cluster::get_id() const {
    return chrom + "_" + std::to_string(s) + "_" + std::to_string(e) ;
}

string Cluster::poa() const {
    uint n_seqs = size();
    vector<string> seqs(n_seqs);
    uint i = 0;
    for (const Fragment &f : fragments) {
        // cout << ">" << i << "\n" << f.seq << endl;
        seqs[i] = f.seq;
        ++i;
    }

    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->out_msa = 0;  // generate Row-Column multiple sequence alignment (RC-MSA)
    abpt->out_cons = 1; // generate consensus sequence
    abpt->disable_seeding = 1;
    // abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
    abpt->progressive_poa = 1;
    abpoa_post_set_para(abpt);

    int *seq_lens = (int*) malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**) malloc(sizeof(uint8_t*) * n_seqs);
    for (uint i = 0; i < n_seqs; ++i) {
        seq_lens[i] = seqs[i].size() ;
        bseqs[i] = (uint8_t*) malloc(sizeof(uint8_t) * seq_lens[i]);
        for (int j = 0; j < seq_lens[i]; ++j) {
            bseqs[i][j] = _nt4_table[(int)seqs[i][j]];
        }
    }

    uint8_t **cons_seq;
    int **cons_cov, *cons_l, cons_n = 0;
    uint8_t **msa_seq;
    int msa_l = 0;

    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);

    // TODO: manage more than 1 consensus
    if (cons_n > 1) {
        cerr << "[W] " << cons_n << " consensus!" << endl;
    }
    // for (int i = 0; i < cons_n; ++i) { }

    string consensus = "";
    for (int j = 0; j < cons_l[0]; ++j) {
        consensus += "ACGTN"[cons_seq[0][j]];
    }

    for (int i = 0; i < cons_n; ++i) {
        free(cons_seq[i]);
        free(cons_cov[i]);
    }
    free(cons_seq);
    free(cons_cov);
    free(cons_l);

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

    return consensus;
}

uint Cluster::get_type() const {
    return fragments.front().t;
}

void Cluster::clear() {
    s = -1;
    e = 0;
    fragments.clear();
}

uint Cluster::size() const {
    return fragments.size();
}

bool Cluster::empty() const {
    return fragments.empty();
}

const Fragment &Cluster::front() const {
    return fragments.front();
}

const Fragment &Cluster::back() const {
    return fragments.back();
}

vector<Fragment>::const_iterator Cluster::begin() const {
    return fragments.begin();
}

vector<Fragment>::const_iterator Cluster::end() const {
    return fragments.end();
}
