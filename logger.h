// Copyright (C) 2024 Maciej Drwal
//
// Permission is granted to copy and distribute verbatim copies and modified
// versions of this file, provided that the copyright notice and this permission
// notice are preserved on all copies and modified versions of this file.
//

#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <string_view>
#include <sstream>
#include <iomanip>

namespace
{
    // Fix for:
    // "The implementation of localtime_s/gmtime_s in Microsoft CRT is incompatible with the C standard
    // since it has reversed parameter order." (https://en.cppreference.com/w/c/chrono/localtime)
    inline void localtime_s_impl(const time_t *t, struct tm *b)
    {
#ifdef _WIN32
        localtime_s(b, t);
#else
        localtime_r(t, b);
#endif
    };
}

namespace utils
{
    enum class Severity
    {
        error,
        warning,
        info,
        debug
    };

    class Logger
    {
    public:
        static Logger &get_instance()
        {
            static Logger instance;
            return instance;
        }

        Logger(const Logger &) = delete;
        void operator=(const Logger &) = delete;
        ~Logger() { m_fout.close(); }

        Logger &set_severity(Severity sev)
        {
            m_severity = sev;
            return *this;
        }

        std::ofstream &operator<<(std::string_view msg)
        {
            m_fout << "\n"
                   << get_tags() << msg;
            if (m_autoflush)
            {
                m_fout << std::flush;
            }
            return m_fout;
        }

    private:
        Logger()
        {
            m_fout.open(m_filename);
            m_fout << get_tags() << "Logger initialized. Log file: " << m_filename;
        }

        std::string_view get_severity() const
        {
            switch (m_severity)
            {
            case Severity::error:
                return "error";
            case Severity::warning:
                return "warning";
            case Severity::debug:
                return "debug";
            case Severity::info:
                return "info";
            default:
                return "unknown severity";
            }
        }

        std::string get_tags() const
        {
            std::stringstream ss;
            const auto t = std::time(nullptr);
            std::tm tm{};
            localtime_s_impl(&t, &tm);
            ss << '[' << std::put_time(&tm, "%F %T%z") << "][" << get_severity() << "] ";
            return ss.str();
        }

        std::string m_filename = "default.log";
        std::ofstream m_fout;
        bool m_autoflush = true;
        Severity m_severity = Severity::info;
    };
}

#define LOG(_s) utils::Logger::get_instance().set_severity(utils::Severity::_s)

#endif
