#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <chrono>

namespace utils
{
    constexpr double ALMOST_ZERO = 1e-7;

    struct Stopwatch
    {
        std::chrono::time_point<std::chrono::steady_clock> t_start;
        void measure_time_start() { t_start = std::chrono::steady_clock::now(); }
        std::string print_elapsed_time()
        {
            auto t_end = std::chrono::steady_clock::now();
            std::stringstream ss;
            ss << "Elapsed time: " << std::chrono::duration_cast<std::chrono::seconds>(t_end - t_start).count() << " s ("
               << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << " ms, "
               << std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count() << " us, "
               << std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count() << " ns)";
            return ss.str();
        }
    };

    template <class T>
    T abs(T x)
    {
        return ((x) < 0.0 ? -(x) : (x));
    }

    template <class T>
    bool is_float_zero(T x)
    {
        return abs(x) < ALMOST_ZERO ? true : false;
    }

    template <class T>
    bool float_equal(T x, T y)
    {
        return is_float_zero(x - y);
    }
}

#define MEASURE_TIME_START(tag) \
    utils::Stopwatch sw_##tag;  \
    sw_##tag.measure_time_start()
#define PRINT_ELAPSED_TIME(tag) sw_##tag.print_elapsed_time()

#endif 
