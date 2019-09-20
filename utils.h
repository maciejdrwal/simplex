#ifndef _UTILS_H_
#define _UTILS_H_

#include <sstream>

#define _abs(x) ((x) < 0.0 ? -(x) : (x))
#define _isfloatzero(x) (_abs(x) < 1e-9 ? true : false)

template<typename T> string tostr(T i)
{
    stringstream ss;
    ss << i;
    return ss.str();
}

#endif 
