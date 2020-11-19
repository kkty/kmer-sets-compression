#ifndef LOG_H_
#define LOG_H_

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

void InitDefaultLogger() {
  spdlog::set_default_logger(spdlog::stderr_color_mt("default"));
}

void EnableDebugLogs() { spdlog::set_level(spdlog::level::debug); }

#endif
