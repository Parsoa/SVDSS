#include <string>
#include <iostream>

#include "config.hpp"

using namespace std ;

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
        ("type", "", cxxopts::value<std::string>())
        ("index", "", cxxopts::value<std::string>())
        ("fastq", "", cxxopts::value<std::string>())
        ("cutoff", "", cxxopts::value<int>())
        ("t,threads", "", cxxopts::value<int>())
        ("a,append", "", cxxopts::value<std::string>())
        ("workdir", "", cxxopts::value<std::string>())
        ("coverage", "", cxxopts::value<int>())
        ("reference", "", cxxopts::value<std::string>())
        ("b,binary", "", cxxopts::value<bool>()->default_value("false"))
        ("aggregate", "", cxxopts::value<bool>())
    ;
}

void Configuration::parse(int argc, char** argv) {
    auto results = parser.parse(argc, argv) ;
    if (results.count("bed")) {
        bed = results["bed"].as<std::string>() ;
    }
    if (results.count("type")) {
        type = results["type"].as<std::string>() ;
    }
    if (results.count("index")) {
        index = results["index"].as<std::string>() ;
    }
    if (results.count("fastq")) {
        fastq = results["fastq"].as<std::string>() ;
    }
    if (results.count("cutoff")) {
        cutoff = results["cutoff"].as<int>() ;
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
    binary = results["binary"].as<bool>() ;
    aggregate = results["aggregate"].as<bool>() ;
    cout << "Coverage: " << coverage << endl ;
}

