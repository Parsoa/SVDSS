#ifndef LPRINT_HPP
#define LPRINT_HPP

#include <initializer_list>
#include <iostream>
#include <string>

using namespace std;

void lprint(
    const initializer_list<string> &tokens, ///< list of strings to print
    int loglevel = 0,  ///< 0 (Info), 1 (Warning), 2 (Error)
    int minlevel = 0); ///< minlevel to print (currently unused)

#endif
