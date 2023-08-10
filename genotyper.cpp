#include "genotyper.hpp"


Genotyper::Genotyper(){
    // according to RAII principle, the constructor should initialize the member variables
    posterior_sv_genotype = {0.0, 0.0, 0.0, 0.0};
}


double Genotyper::prob_sfs_sv(int sfs_length, int sfs, allele_t a){

    (void)(sfs_length); // suppress unused parameter warning

    if (sfs == 1 && a == 1){
        return 0.8;
    }
    else if (sfs == 0 && a == 0){
        return 0.95;
    }
    else if (sfs == 0 && a == 1){
        return 0.2;
    }
    else{ // sfs == 1 && a == 0
        return 0.05;
    }

#ifdef DEBUG
    assert(false); // should never reach here
#endif

    return -1;
}


double Genotyper::priorHap(int hap, int h){

    if (h == 3){
        return 0.5;
    }

    if (hap == h){
        return 0.95;
    }
    else{
        return 0.05;
    }

#ifdef DEBUG
    assert(false); // should never reach here
#endif

    return -1;
}


double Genotyper::prior_genotype(variant_t v){

    if (v == make_variant(0, 0)){
        return 0.5;
    }
    else if (v == make_variant(0, 1) || v == make_variant(1, 0)){
        return 0.22;
    }
    else{ // v == make_variant(1, 1)
        return 0.06;
    }

#ifdef DEBUG
    assert(false); // should never reach here
#endif
    
    return -1;
}


double Genotyper::likelihood_read_give_gen_hap(int s, int hap, allele_t a){

    (void)(hap); // suppress unused parameter warning

    return this->prob_sfs_sv(20, s, a);
}


double Genotyper::likelihood_read_give_genotype(read_t r, variant_t v){

    int hap1 = 1;
    int hap2 = 2;

    double read_likelihood_allel_0 = this->likelihood_read_give_gen_hap(std::get<0>(r), hap1, std::get<0>(v)) * this->priorHap(hap1, std::get<1>(r));
    double read_likelihood_allel_1 = this->likelihood_read_give_gen_hap(std::get<0>(r), hap2, std::get<1>(v)) * this->priorHap(hap2, std::get<1>(r));

    return read_likelihood_allel_0 + read_likelihood_allel_1;
}


double Genotyper::likelihood_allreads_give_genotype(std::vector<read_t> r_vec, variant_t v ){

    double log_posterior_prob = 0.0;

    for (auto r : r_vec){
        log_posterior_prob += this->likelihood_read_give_genotype(r, v);
    }

    return log_posterior_prob;
}


/**
 * @brief compute the posterior probability of each genotype given all reads
 * 
 * @param r_vec a vector of reads.
 *        Each read is a tuple of two integers:
 *        - the first integer is indication whether the read has a SFS supporting the variant allele (1) or not (0)
 *        - the second integer is the predicted allele (haplo?)tag of the read (1 for haplotype 1, 2 for haplotype 2, 3 for NOTAG)    
 * @return int 0 if success
 */
int Genotyper::posterior_sv_genotype_give_reads(std::vector<read_t> r_vec){

    int i=0;
    int total = 0;

    std::vector<std::tuple<int,int> > g = {std::make_tuple(0, 0), std::make_tuple(0, 1), std::make_tuple(1, 0), std::make_tuple(1, 1)};

    for (auto v : g){
        this->posterior_sv_genotype[i] = this->likelihood_allreads_give_genotype(r_vec, v) + log(this->prior_genotype(v));
        total += exp(this->posterior_sv_genotype[i]);
        i++;
    }

    for (i=0; i<4; i++){
        this->posterior_sv_genotype[i] = exp(this->posterior_sv_genotype[i]) / total;
    }

    return 0;
}


std::vector<double> Genotyper::get_posterior_sv_genotype(){

    return this->posterior_sv_genotype;
}