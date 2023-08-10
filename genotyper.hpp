#ifndef GENOTYPER_HPP
#define GENOTYPER_HPP

#include <string>
#include <tuple>
#include <functional>
#include <vector>
#include <cmath>


typedef std::tuple<int, int> variant;
typedef std::tuple<int, int> read;
typedef int allele;

static variant make_variant(int a, int b){return std::make_tuple(a, b);}
static read make_read(int a, int b){return std::make_tuple(a, b);}

// Define a class "Genotyper" that will be used to compute the genotype of Structural Variants
// the class has the following private member functions:
// 1. prob_sfs_sv: compute the likelihood to observe a read based on it has a SFS and given the genotype of the SV that read spans
// 2. priorHap: compute the prior probability to observe a read based on its predicted haplotype
// 3. prior_genotype: define prior probability of each genotype
// 4. likelihood_read_give_gen_hap
// 5. likelihood_allreads_give_genotype: log likelihood of all reads
// this class has the following public member functions:
// 1. posterier_sv_genotype_give_reads: compute the posterior probability of each genotype given all reads
// 2. get_posterior_sv_genotype: compute the posterior probability of each genotype
// this class has the following private member variables:
// 1. vector of posterior probability of each genotype
class Genotyper{
private:
    //
    double prob_sfs_sv(int sfs_length, int sfs, allele a);
    // 
    double priorHap(int hap, int h);
    // 
    double prior_genotype(variant v);
    //
    double likelihood_read_give_gen_hap(int s, int hap, allele a);
    // 
    double likelihood_read_give_genotype(read r, variant v);
    //
    double likelihood_allreads_give_genotype(std::vector<read> r_vec, variant v );

public:
    Genotyper();
    //
    int posterior_sv_genotype_give_reads(std::vector<read> r_vec);
    //
    std::vector<double> get_posterior_sv_genotype();

private:
    std::vector<double> posterior_sv_genotype;
};



#endif // GENOTYPER_HPP