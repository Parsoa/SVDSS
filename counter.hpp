#ifndef CNT_HPP
#define CNT_HPP

#include <string>
#include <unordered_map>

class Counter {

public:

    void load_counts(std::string path) ;

    std::unordered_map<std::string, int> counts ;

} ;

#endif
