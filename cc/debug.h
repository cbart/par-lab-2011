// Debug tools.
// author: Cezary Bartoszuk <cbart@students.mimuw.edu.pl>
//     id: cb277617@students.mimuw.edu.pl

#ifndef _DEBUG_TOOLS_H_
#define _DEBUG_TOOLS_H_

#include <iostream>

namespace debug
{

    typedef unsigned long debug_level_t;

    #ifdef DEBUGLEVEL
        const debug_level_t DEBUG_LEVEL = DEBUGLEVEL;
    #else
        const debug_level_t DEBUG_LEVEL = 0;
    #endif

    const std::string DIAG_PROG_NAME = "DEBUG";

    #ifdef ERRORLEVEL
        const debug_level_t ERROR_DEBUG_LEVEL = ERRORLEVEL;
    #else
        const debug_level_t ERROR_DEBUG_LEVEL = 2;
    #endif

    #ifdef WARNINGLEVEL
        const debug_level_t WARNING_DEBUG_LEVEL = WARNINGLEVEL;
    #else
        const debug_level_t WARNING_DEBUG_LEVEL = 2;
    #endif

    #ifdef INFOLEVEL
        const debug_level_t INFO_DEBUG_LEVEL = INFOLEVEL;
    #else
        const debug_level_t INFO_DEBUG_LEVEL = 2;
    #endif

    // A tool for managing diagnostic info.
    class diagstream
    {

        private:

            // Minimal level of debugging with which
            // this object will print the output.
            const debug_level_t MIN_DEBUG_LEVEL;

        public:

            typedef std::ostream& (*message_fun_t)(std::ostream&);

            // Creates new diagnostic stream with defined
            // `min_debug_lv` (minimal debug level).
            diagstream(debug_level_t min_debug_lv)
                : MIN_DEBUG_LEVEL(min_debug_lv)
            {
            }

            // 'Send to stream' operator designed for all kinds
            // of messages.
            template<typename T>
                friend diagstream& operator<<(diagstream& ds, const T& message);

            // 'Send to stream' operator designed for operations
            // like std::flush, std::endl, etc...
            friend diagstream& operator<<
                (diagstream& ds, message_fun_t message_fun);

    };

    template<typename T>
    inline static diagstream& operator<<(diagstream& ds, const T& message)
    {
        if(DEBUG_LEVEL >= ds.MIN_DEBUG_LEVEL)
            ::std::cerr << message;
        return ds;
    }

    inline diagstream& operator<<
    (diagstream& ds, diagstream::message_fun_t message_fun)
    {
        if(DEBUG_LEVEL >= ds.MIN_DEBUG_LEVEL)
            ::std::cerr << message_fun;
        return ds;
    }

    // Returns stream for error logging.
    inline diagstream& get_errlog_stream(const std::string& info_str)
    {
        static diagstream error_stream = diagstream(ERROR_DEBUG_LEVEL);
        if(DEBUG_LEVEL >= ERROR_DEBUG_LEVEL)
            ::std::cerr << info_str;
        return error_stream;
    }

    // Returns stream for warnings logging.
    inline diagstream& get_warnlog_stream(const std::string& info_str)
    {
        static diagstream warning_stream = diagstream(WARNING_DEBUG_LEVEL);
        if(DEBUG_LEVEL >= WARNING_DEBUG_LEVEL)
            ::std::cerr << info_str;
        return warning_stream;
    }

    // Returns stream for information logging.
    inline diagstream& get_infolog_stream(const std::string& info_str)
    {
        static diagstream info_stream = diagstream(INFO_DEBUG_LEVEL);
        if(DEBUG_LEVEL >= INFO_DEBUG_LEVEL)
            ::std::clog << info_str;
        return info_stream;
    }

    // Some useful defines (so we won't always write __PRETTY_FUNCTION__ or
    // "(II) Quaternion::Quaternion(const Quaternion& q); "

    // Indicator strings.
    static const std::string DBG_ERROR_INDICATOR("(EE) ");
    static const std::string DBG_WARNING_INDICATOR("(WW) ");
    static const std::string DBG_INFO_INDICATOR("(II) ");
    static const std::string DBG_STRING_WHITESPACE(" ");

    // `err` definition.
    #define err get_errlog_stream(debug::DBG_ERROR_INDICATOR \
            + __PRETTY_FUNCTION__ + debug::DBG_STRING_WHITESPACE)

    // `warn` definition.
    #define warn get_warnlog_stream(debug::DBG_WARNING_INDICATOR \
            + __PRETTY_FUNCTION__ + debug::DBG_STRING_WHITESPACE)

    // `info` definition.
    #define info get_infolog_stream(debug::DBG_INFO_INDICATOR \
            + __PRETTY_FUNCTION__ + debug::DBG_STRING_WHITESPACE)

    void print_self_sig(const std::string& sig )
    {
        ::std::clog << DBG_INFO_INDICATOR << sig << std::endl;
    }

    // Prints function signature to std::clog;
    #define print_sig print_self_sig(__PRETTY_FUNCTION__)

}

#endif
