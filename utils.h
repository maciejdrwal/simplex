#ifndef UTILS_H
#define UTILS_H

#include <sstream>

namespace utils
{
    template <class T>
    T abs(T x)
    {
        return ((x) < 0.0 ? -(x) : (x));
    }

    template <class T>
    bool is_float_zero(T x)
    {
        return abs(x) < 1e-9 ? true : false;
    }

    template <class T>
    std::string to_str(T i)
    {
        std::stringstream ss;
        ss << i;
        return ss.str();
    }
}

#endif 
