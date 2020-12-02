#ifndef LOG_H_
#define LOG_H_

#include <memory>

#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

void InitDefaultLogger() {
  std::shared_ptr<spdlog::logger> logger = spdlog::stderr_color_mt("default");
  logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] [%t] %v");
  spdlog::set_default_logger(logger);
}

void EnableDebugLogs() { spdlog::set_level(spdlog::level::debug); }

#endif
