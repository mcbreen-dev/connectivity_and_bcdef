/*─────────────────────────────────────────────────────────────
  File: src/logger.cpp

  Minimal console + file logger.

  This module implements the Logger class used throughout the project
  to emit user-facing status messages, warnings, and errors.

  Features:
    - Writes to stdout (always) and optionally to a log file
    - Thread-safe: write() is protected by a mutex
    - ANSI color tags for interactive terminals only
        * INFO  -> green background
        * WARN  -> yellow background
        * ERROR -> red background
      When stdout is not a TTY (e.g., redirected to a file), the logger
      falls back to plain, uncolored tags so output remains readable.

  Notes:
    - Color detection uses isatty(fileno(stdout)) from unistd.h.
    - File logs always use plain tags (no escape codes).
─────────────────────────────────────────────────────────────*/
#include "logger.hpp"
#include <unistd.h>          // isatty

namespace bcdef {

/*=====================================================================
  ANSI formatting constants

  These are used only for terminal output (stdout) when stdout is a TTY.
  They are *not* used for file output, to avoid polluting log files with
  escape codes.
=====================================================================*/
#define BG_GRN  "\033[102m"
#define BG_YEL  "\033[103m"
#define BG_RED  "\033[101m"
#define RESET   "\033[0m"

/*=====================================================================
  bare_tag

  Return a plain (uncolored) log level tag for a Logger::Level.

  This is used:
    - for file logging (always plain)
    - for console logging when stdout is not a terminal (isatty == false)

  Output strings are short and stable to make grepping logs easy.
=====================================================================*/
static const char* bare_tag(Logger::Level l)      // plain for file log
{
    switch (l) {
        case Logger::Level::INFO: return "INFO";
        case Logger::Level::WARN: return "WARN";
        case Logger::Level::ERR : return "ERROR";
    }
    return "UNKWN";
}

/*=====================================================================
  Logger::level_tag  (static)

  Return the level tag string that should be used for console output.

  Behavior:
    - If stdout is not a terminal (not interactive), return a plain tag
      (no ANSI escapes).
    - Otherwise, return a colored tag using background color escapes.

  Important:
    - The returned pointer refers to string literals (static storage),
      so it remains valid for the lifetime of the program.
=====================================================================*/
const char* Logger::level_tag(Level l)
{
    /* color only if writing to an interactive terminal         */
    if (!isatty(fileno(stdout))) return bare_tag(l);

    switch (l) {
        case Level::INFO: return BG_GRN "INFO" RESET;
        case Level::WARN: return BG_YEL "WARN" RESET;
        case Level::ERR : return BG_RED "ERROR" RESET;
    }
    return "UNKWN";
}

/*=====================================================================
  Logger lifecycle

  The logger optionally writes to a file:
    - Constructor opens the file in truncation mode (fresh log per run)
    - Destructor closes the file if it was opened successfully

  Console output is always written, regardless of whether the file opens.
=====================================================================*/
Logger::Logger(const std::string& p){ file_.open(p,std::ios::trunc); }
Logger::~Logger(){ if(file_) file_.close(); }

/*=====================================================================
  Logger::write

  Emit one log line to:
    - stdout (colored tag if interactive terminal)
    - file (if open), using plain tags only

  Thread safety:
    - A mutex protects both stdout and file writes so that concurrent
      writers do not interleave partial lines.

  Formatting:
    "<TAG>  <message>\n"
=====================================================================*/
void Logger::write(Level lvl,const std::string& msg)
{
    std::string line = std::string(level_tag(lvl)) + "  " + msg + '\n';
    std::lock_guard<std::mutex> g(mtx_);
    std::cout << line << std::flush;         // colorised
    if (file_) file_ << bare_tag(lvl) << "  " << msg << '\n';
}

} // namespace bcdef
