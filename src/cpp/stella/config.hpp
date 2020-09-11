#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>
#include <vector>

#include "cxxopts.hpp"
#include "bed_utils.hpp"

class Configuration {

private:
    static Configuration* instance ;

public:
    static Configuration* getInstance() ;

    void parse(int argc, char* argv[]) ;

    int threads ;
    int batch_size = 1000 ;

    std::string bed ;
    std::string workdir ;
    std::string reference ;

private:

    Configuration() ;

    Configuration(Configuration const&) = delete ;
    void operator=(Configuration const&) = delete ;

    Configuration& operator[](std::string) ;
    
    cxxopts::Options parser ;
};

#endif
