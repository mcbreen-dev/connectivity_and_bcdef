/*─────────────────────────────────────────────────────────────
  File: src/logger.cpp
─────────────────────────────────────────────────────────────*/
#include "logger.hpp"
#include <unistd.h>          // isatty

namespace fs {

/* ANSI background colours */
#define BG_GRN  "\033[102m"
#define BG_YEL  "\033[103m"
#define BG_RED  "\033[101m"
#define RESET   "\033[0m"

static const char* bare_tag(Logger::Level l)      // plain for file log
{
    switch (l) {
        case Logger::Level::INFO: return "INFO";
        case Logger::Level::WARN: return "WARN";
        case Logger::Level::ERR : return "ERROR";
    }
    return "UNKWN";
}

/* ---------- THE REQUIRED CLASS STATIC ---------------------- */
const char* Logger::level_tag(Level l)
{
    /* colour only if writing to an interactive terminal         */
    if (!isatty(fileno(stdout))) return bare_tag(l);

    switch (l) {
        case Level::INFO: return BG_GRN "INFO" RESET;
        case Level::WARN: return BG_YEL "WARN" RESET;
        case Level::ERR : return BG_RED "ERROR" RESET;
    }
    return "UNKWN";
}

/* unchanged ctor / dtor / write() ---------------------------- */
Logger::Logger(const std::string& p){ file_.open(p,std::ios::trunc); }
Logger::~Logger(){ if(file_) file_.close(); }

void Logger::write(Level lvl,const std::string& msg)
{
    std::string line = std::string(level_tag(lvl)) + "  " + msg + '\n';
    std::lock_guard<std::mutex> g(mtx_);
    std::cout << line << std::flush;         // colourised
    if (file_) file_ << bare_tag(lvl) << "  " << msg << '\n';
}

} // namespace fs
