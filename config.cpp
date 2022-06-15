#include "config.hpp"

Configuration* Configuration::instance = nullptr ;

Configuration* Configuration::getInstance() {
    if (instance == nullptr) {
        instance = new Configuration() ;
    }
    return instance ; 
}

Configuration::Configuration() :
    parser("Stella, mapping-free variation discovery.") {
    parser.add_options()
        ("bed", "", cxxopts::value<std::string>())
        ("bam", "", cxxopts::value<std::string>())
        ("vcf", "", cxxopts::value<std::string>())
        ("type", "", cxxopts::value<std::string>())
        ("index", "", cxxopts::value<std::string>())
        ("fasta", "", cxxopts::value<std::string>())
        ("fastq", "", cxxopts::value<std::string>())
        ("target", "", cxxopts::value<std::string>())
        ("sfsbam", "", cxxopts::value<std::string>())
        ("append", "", cxxopts::value<std::string>())
        ("workdir", "", cxxopts::value<std::string>())
        ("reference", "", cxxopts::value<std::string>())
        ("cutoff", "", cxxopts::value<int>())
        ("batches", "", cxxopts::value<int>())
        ("overlap", "", cxxopts::value<int>())
        ("coverage", "", cxxopts::value<int>())
        ("t,threads", "", cxxopts::value<int>())
        ("min-sv-length", "", cxxopts::value<int>())
        ("min-cluster-weight", "", cxxopts::value<int>())
        ("clipped", "", cxxopts::value<bool>()->default_value("false"))
        ("assemble", "", cxxopts::value<bool>()->default_value("false"))
        ("putative", "", cxxopts::value<bool>()->default_value("true"))
        ("b,binary", "", cxxopts::value<bool>()->default_value("false"))
        ("aggregate", "", cxxopts::value<bool>()->default_value("false"))
        ("selective", "", cxxopts::value<bool>()->default_value("true"))
        ("version", "Print version information.")
        ("help", "Print this help.")
    ;
}

void Configuration::parse(int argc, char** argv) {
    auto results = parser.parse(argc, argv) ;
    if (results.count("vcf")) {
        vcf = results["vcf"].as<std::string>() ;
    }
    bam = "" ;
    if (results.count("bam")) {
        bam = results["bam"].as<std::string>() ;
    }
    if (results.count("sfsbam")) {
        sfsbam = results["sfsbam"].as<std::string>() ;
    }
    if (results.count("bed")) {
        bed = results["bed"].as<std::string>() ;
    }
    if (results.count("type")) {
        type = results["type"].as<std::string>() ;
    }
    if (results.count("index")) {
        index = results["index"].as<std::string>() ;
    }
    fastq = "" ;
    if (results.count("fastq")) {
        fastq = results["fastq"].as<std::string>() ;
    }
    fasta = "" ;
    if (results.count("fasta")) {
        fasta = results["fasta"].as<std::string>() ;
    }
    target = "" ;
    if (results.count("target")) {
        target = results["target"].as<std::string>() ;
    }
    if (results.count("cutoff")) {
        cutoff = results["cutoff"].as<int>() ;
    }
    if (results.count("overlap")) {
        overlap = results["overlap"].as<int>() ;
    }
    if (results.count("threads")) {
        threads = results["threads"].as<int>() ;
    }
    if (results.count("workdir")) {
        workdir = results["workdir"].as<std::string>() ;
    } else {
        workdir = "." ;
    }
    if (results.count("append")) {
        append = results["append"].as<std::string>() ;
    } else {
        append = "" ;
    }
    if (results.count("coverage")) {
        coverage = results["coverage"].as<int>() ;
    }
    if (results.count("reference")) {
        reference = results["reference"].as<std::string>() ;
    }
    if (results.count("batches")) {
        aggregate_batches = results["batches"].as<int>() ;
    }
    if (results.count("min-sv-length")) {
        min_sv_length = max(25, results["min-sv-length"].as<int>()) ;
    }
    if (results.count("min-cluster-weight")) {
        min_cluster_weight = results["min-cluster-weight"].as<int>() ;
    }
    binary = results["binary"].as<bool>() ;
    clipped = results["clipped"].as<bool>() ;
    assemble = results["assemble"].as<bool>() ;
    putative = results["putative"].as<bool>() ;
    aggregate = results["aggregate"].as<bool>() ;
    selective = results["selective"].as<bool>() ;
    version = results["version"].as<bool>();
    help = results["help"].as<bool>();
}
