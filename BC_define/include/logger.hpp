/*─────────────────────────────────────────────────────────────
  File: include/logger.hpp
─────────────────────────────────────────────────────────────*/
#pragma once

#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <chrono>
#include <iomanip>

namespace fs {

class Logger
{
public:
    enum class Level { INFO, WARN, ERR };

    explicit Logger(const std::string& logfile_path);
    ~Logger();

    /** Thread-safe write to log file **and** stdout */
    void write(Level level, const std::string& msg);

    /* convenience wrappers */
    void info (const std::string& m) { write(Level::INFO , m); }
    void warn (const std::string& m) { write(Level::WARN , m); }
    void error(const std::string& m) { write(Level::ERR  , m); }

    static const char* level_tag(Level);

private:
    std::ofstream file_;
    std::mutex    mtx_;

    static std::string timestamp();
};

} // namespace fs