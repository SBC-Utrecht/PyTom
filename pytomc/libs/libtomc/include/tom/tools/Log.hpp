/***************************************************************************//**
 * \file Log.hpp
 * \brief Class for logging data.
 * \author  Michael Stoelken
 * \version 0.1
 * \date    29.01.2009
 *
 * Class for logging data.
 ******************************************************************************/
#ifndef ___INCLUDE__TOM__TOOLS__LOG_HPP__
#define ___INCLUDE__TOM__TOOLS__LOG_HPP__




#include <ostream>
#include <fstream>
#include <iomanip>

#include <cerrno>
#include <cassert>
#include <sstream>
#include <stdexcept>


#include <tom/tools/ostream_swapper.hpp>
#include <tom/tools/snippets.hpp>



/// Use this define to wrapp calls of the logging which are done with an
/// a priory known log-level. This can be used to cancel the entire call of
/// the logging, if the log-level is too low.\n
/// This is useful for statements like 'logger(level, tag) << some_expensive_call();'\n
/// Beware of side effects!!
#define TLOG(level, cmd)                                                        \
    do {                                                                        \
        if ( (level) <= tom::tools::Log::max_log_level) {                       \
            cmd ;                                                               \
        }                                                                       \
    } while (0)

#define TLOG_ERROR(cmd)             TLOG(tom::tools::Log::ERROR  , cmd)
#define TLOG_WARNING(cmd)           TLOG(tom::tools::Log::WARNING, cmd)
#define TLOG_INFO(cmd)              TLOG(tom::tools::Log::INFO   , cmd)
#define TLOG_DEBUG(cmd)             TLOG(tom::tools::Log::DEBUG  , cmd)
#define TLOG_DEBUGV(cmd)            TLOG(tom::tools::Log::DEBUGV , cmd)

namespace tom {
namespace tools {




struct s_Endl {
    s_Endl() { }
};
extern const s_Endl endl;


/***************************************************************************//**
 * Logging class.
 ******************************************************************************/
class Log {

public:
    enum Level {  ERROR
                , WARNING
                , INFO
                , DEBUG
                , DEBUGV
                };

    Log(std::ostream *s_default, char separator=';');
    ~Log();

    void redirect(const std::string &filename);
    void redirect(std::ostream &stream);
    template<typename T> void log(Level level, const std::string &msgid, const T &s);
    Log &log(Level level, const std::string &msgid);
    template<typename T> Log &operator<<(T e);
    static std::string Level2str(Level level);
    void newline();

    void setTagWidth(int setw_tag);

    void flush();

    #if defined(TOM_TOOLS_MAX_LOG_LEVEL)
    static const Log::Level max_log_level = TOM_TOOLS_MAX_LOG_LEVEL;
    #else
    #if defined(NDEBUG)
    static const Log::Level max_log_level = Log::INFO;
    #else
    static const Log::Level max_log_level = Log::DEBUGV;
    #endif
    #endif
private:
    Log(const Log &): s_default_(0), s_(0), swp_(0)    { assert(!"HIDDEN"); }
    Log &operator=(const Log &)         { assert(!"HIDDEN"); return *this; }

    void make_new_entry(Level level, const std::string &msgid);
    void make_new_entry();

    std::ostream *s_default_;
    std::ostream *s_;
    tom::tools::ostream_swapper *swp_;
    bool entry_closed_;
    bool is_redirected_;
    Level old_level_;
    std::string old_msgid_;
    int setw_tag_;
    char separator_;
};


} // namespace tools
} // namespace tom









/***************************************************************************//**
 *
 ******************************************************************************/
template<typename T>
inline void tom::tools::Log::log(Level level, const std::string &msgid, const T &s) {
    make_new_entry(level, msgid);
    if (old_level_ <= max_log_level) {
        *s_ << s << '\n';
        entry_closed_ = true;
    }
}




/***************************************************************************//**
 *
 ******************************************************************************/
inline tom::tools::Log &tom::tools::Log::log(Level level, const std::string &msgid) {
    make_new_entry(level, msgid);
    return *this;
}



namespace tom {
namespace tools {
/***************************************************************************//**
 *
 ******************************************************************************/
template<typename T> Log &Log::operator<<(T e) {
    if (old_level_ <= max_log_level) {
        if (entry_closed_) {
            make_new_entry();
        }
        *s_ << e;
    }
    return *this;
}
/***************************************************************************//**
 *
 ******************************************************************************/
template<>
inline Log &Log::operator<<(s_Endl) {
    if (old_level_ <= max_log_level) {
        entry_closed_ = true;
        *s_ << '\n';
    }
    return *this;
}
} // namespace tools
} // namespace tom


/***************************************************************************//**
 *
 ******************************************************************************/
inline void tom::tools::Log::newline() {
    if (old_level_ <= max_log_level) {
        *s_ << '\n';
        entry_closed_ = true;
    }
}


/***************************************************************************//**
 *
 ******************************************************************************/
inline void tom::tools::Log::flush() {
    (*s_).flush();
}



/***************************************************************************//**
 *
 ******************************************************************************/
inline std::string tom::tools::Log::Level2str(Level level) {
    switch (level) {
        case INFO:
            return "INFO";
        case ERROR:
            return "ERROR";
        case WARNING:
            return "WARNING";
        case DEBUG:
            return "DEBUG";
        case DEBUGV:
            return "DEBUGV";
    }
    assert(!"INVALID LOG LEVEL");
    return "";
}




/***************************************************************************//**
 *
 ******************************************************************************/
inline void tom::tools::Log::make_new_entry(Level level, const std::string &msgid) {
    old_level_ = level;
    old_msgid_ = msgid;
    make_new_entry();
}



/***************************************************************************//**
 *
 ******************************************************************************/
inline void tom::tools::Log::setTagWidth(int setw_tag) {
    setw_tag_ = setw_tag;
}



/***************************************************************************//**
 *
 ******************************************************************************/
inline void tom::tools::Log::make_new_entry() {
    if (old_level_ <= max_log_level) {
        if (!entry_closed_) {
            *s_ << '\n';
        } else {
            entry_closed_ = false;
        }
        *s_            << std::setw(       19) << tom::tools::now2str() << ' ' << separator_ << ' ' <<
            std::right << std::setw(        7) << Level2str(old_level_) << ' ' << separator_ << ' ' <<
            std::left  << std::setw(setw_tag_) << old_msgid_            << ' ' << separator_ << ' ';
    }
}



#endif

