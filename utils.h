#ifndef _UTILS_H_
#define _UTILS_H_

#include <sstream>

namespace utils
{
    template <class T>
    T abs(T x)
    {
        return ((x) < 0.0 ? -(x) : (x));
    }

    template <class T>
    bool isfloatzero(T x)
    {
        return abs(x) < 1e-9 ? true : false;
    }

    template <class T>
    std::string tostr(T i)
    {
        std::stringstream ss;
        ss << i;
        return ss.str();
    }
}

#endif 
