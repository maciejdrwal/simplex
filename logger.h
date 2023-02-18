#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <fstream>
#include <string_view>
#include <sstream>
#include <iomanip>

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

        static Logger & get_instance()
        {
            static Logger instance;
            return instance;
        }

        Logger(const Logger &) = delete;
        void operator=(const Logger &) = delete;
        ~Logger() { m_fout.close(); }

        Logger & set_severity(Severity sev)
        {
            m_severity = sev;
            return *this;
        }

private:
        
        Logger() 
        { 
            m_fout.open(m_filename); 
            m_fout << get_tags() << "Logger initialized. Log file: " << m_filename;
        }
        
        std::string_view get_severity() const
        {
            switch(m_severity)
            {
                case Severity::error: return "error";
                case Severity::warning: return "warning";
                case Severity::debug: return "debug";
                default: return "unknown severity";
            }
        }

        std::string get_tags() const
        {
            std::stringstream ss;
            const auto t = std::time(nullptr);
            ss << '[' << std::put_time(std::localtime(&t), "%F %T%z") << "][" << get_severity() << "] ";
            return ss.str();
        }

        std::ofstream & operator<<(const std::string & msg)
        {
            m_fout << "\n" << get_tags() << msg;
            if (m_autoflush)
            {
                m_fout << std::flush;
            }
            return m_fout;
        }

    private:

        std::string m_filename = "default.log";
        std::ofstream m_fout;
        bool m_autoflush = true;
        Severity m_severity = Severity::info;
    };
}

#define LOG(_s) utils::Logger::get_instance().set_severity(utils::Severity::_s)

#endif
