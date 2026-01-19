/*─────────────────────────────────────────────────────────────
  File: include/logger.hpp

  Logger: small thread-safe logging utility.

  This class is used in:
    - bc_define_main.cpp (CLI orchestration + fatal error reporting)
    - connectivity_core.cpp (progress + connectivity/boundary reporting)
    - mesh_utils.cpp (purge/prune reporting)
    - bc_define_cgns.cpp (IO benchmark + BC write stage reporting)
    - bc_dump.cpp (open/finish/error reporting)
    - resolve_matches_* code paths (via caller logs)

  Output behavior:
    - Writes a formatted line to stdout
    - Writes a plain (non-ANSI) line to a log file
    - Guards both outputs with a mutex to prevent interleaving

  Notes about sizes:
    - std::mutex is typically ~40 bytes on many libstdc++ builds (platform dependent)
    - std::ofstream is a non-trivial object that owns a file descriptor/buffer
    - Logger is intended to be a process-level utility object (constructed once per tool)
─────────────────────────────────────────────────────────────*/
#pragma once

#include <fstream>    // std::ofstream (file output stream)
#include <iostream>   // std::cout (console output stream)
#include <mutex>      // std::mutex, std::lock_guard (thread synchronization)
#include <string>     // std::string (message storage/formatting)

namespace fs {

class Logger
{
public:
    /*------------------------------------------------------------------
      Level: severity enum for log messages.

      Underlying type:
        - not explicitly specified; defaults to int (typically 4 bytes)

      Values:
        - INFO : normal status/progress messages
        - WARN : recoverable issues or noteworthy conditions
        - ERR  : errors; used in catch blocks and fatal paths
    ------------------------------------------------------------------*/
    enum class Level { INFO, WARN, ERR };

    /*------------------------------------------------------------------
      Construct a Logger that writes to a file path.

      Parameter:
        logfile_path : filesystem path (UTF-8 string)

      Effects:
        - Opens/truncates the log file (implementation uses std::ios::trunc)
        - After construction, write() is valid

      Lifetime:
        - The Logger object is stack-allocated in main() for bc_define and bc_dump
        - Destructor closes the file if open
    ------------------------------------------------------------------*/
    explicit Logger(const std::string& logfile_path);
    ~Logger(); // Destructor. Closes the file stream if open.

    /*------------------------------------------------------------------
      write: thread-safe logging to both stdout and the log file.

      Parameters:
        level : Logger::Level
        msg   : message text (already formatted by caller; no printf-style args)

      Output formatting (see src/logger.cpp):
        - stdout: may include ANSI color codes if stdout is a TTY
        - file  : always plain tags (INFO/WARN/ERROR), no ANSI sequences

      Thread-safety:
        - Uses std::lock_guard<std::mutex> to serialize output
        - Prevents line interleaving when called from multiple threads
          (connectivity matching can run multithreaded; logging remains readable)
    ------------------------------------------------------------------*/
    void write(Level level, const std::string& msg);

    /*------------------------------------------------------------------
      Convenience wrappers for common severities.

      These forward to write(...) with fixed Level.
      They do not add formatting beyond severity selection.
    ------------------------------------------------------------------*/
    void info (const std::string& m) { write(Level::INFO , m); }
    void warn (const std::string& m) { write(Level::WARN , m); }
    void error(const std::string& m) { write(Level::ERR  , m); }

    /*------------------------------------------------------------------
      level_tag: maps severity to a label string.

      Return:
        const char* pointing to a static string.

      Implementation details (src/logger.cpp):
        - If stdout is not a TTY (isatty(stdout)==false):
            returns plain "INFO"/"WARN"/"ERROR"
        - If stdout is a TTY:
            returns strings containing ANSI background color codes

      Caller usage:
        - write() uses level_tag() when formatting the stdout line
        - File output uses a plain (non-ANSI) tag function in logger.cpp
    ------------------------------------------------------------------*/
    static const char* level_tag(Level);

private:
    /*------------------------------------------------------------------
      file_: output file stream.

      Type:
        std::ofstream

      Ownership:
        - Owns the underlying file descriptor/buffer when open
        - Opened in the constructor
        - Closed in the destructor (if open)

      Note:
        - std::ofstream is not thread-safe by itself; access is guarded by mtx_.
    ------------------------------------------------------------------*/
    std::ofstream file_;

    /*------------------------------------------------------------------
      mtx_: mutex guarding stdout and file writes.

      Type:
        std::mutex

      Role:
        - Ensures each log line is emitted as an atomic unit
          (single contiguous line to stdout, and single line to file)
    ------------------------------------------------------------------*/
    std::mutex    mtx_;
};

} // namespace fs