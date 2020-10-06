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
        ("threads", "", cxxopts::value<int>())
        ("workdir", "", cxxopts::value<std::string>())
        ("reference", "", cxxopts::value<std::string>())
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
    if (results.count("threads")) {
        threads = results["threads"].as<int>() ;
    }
    if (results.count("workdir")) {
        workdir = results["workdir"].as<std::string>() ;
    }
    if (results.count("reference")) {
        reference = results["reference"].as<std::string>() ;
    }
}

