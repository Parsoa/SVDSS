#include "lprint.hpp"

void lprint(const initializer_list<string> &tokens, int loglevel, int minlevel) {
    if (loglevel < minlevel) {
        return;
    }
    char ctype;
    switch (loglevel) {
        case 0:
            ctype = 'I';
            break;
        case 1:
            ctype = 'W';
            break;
        case 2:
            ctype = 'E';
            break;
        default:
            ctype = 'I';
    }
    cerr << "[" << ctype << "]" << " " ;
    for (const string &token : tokens) {
        cerr << token << " ";
    }
    cerr << endl;
}
