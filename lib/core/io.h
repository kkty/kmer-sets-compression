#ifndef CORE_IO_H_
#define CORE_IO_H_

#include <cstdio>
#include <fstream>
#include <istream>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"

// Opens a file and reads lines in it. If "decompressor" is not empty, it will
// be used to decompress the file.
// Example:
//   ... = ReadLiens("foo.txt", "");
//   ... = ReadLines("foo.txt.bz2", "bzip2 -d");
absl::StatusOr<std::vector<std::string>> ReadLines(
    const std::string& file_name, const std::string& decompressor) {
  if (decompressor.empty()) {
    std::ifstream file(file_name);

    if (file.fail()) {
      return absl::InternalError("failed to open file");
    }

    std::vector<std::string> lines;

    std::string s;
    while (std::getline(file, s)) {
      lines.push_back(s);
    }

    return lines;
  }

  std::FILE* f =
      popen(absl::StrFormat("%s < %s", decompressor, file_name).c_str(), "r");

  if (f == NULL) {
    return absl::InternalError("failed to open a sub-process");
  }

  std::string s;

  // Reads data from the process.
  {
    const int n = 8192;
    char buf[n];

    while (fgets(buf, n, f) != NULL) {
      s += std::string(buf);
    }
  }

  // Waits for the process to exit.
  {
    const int exit_status = pclose(f);

    if (exit_status != 0) {
      return absl::InternalError(absl::StrFormat(
          "process failed with non-zero exit code: %d", exit_status));
    }
  }

  if (!s.empty() && s.back() == '\n') s.pop_back();

  std::vector<std::string> lines = absl::StrSplit(s, '\n');

  return lines;
}

// Writes lines to a file. If "compressor" is not empty, it will be used to
// compress the file.
// Example:
//   ... = WriteLines("foo.txt", "", lines);
//   ... = WriteLines("foo.txt.bz2", "bzip2", lines);
absl::Status WriteLines(const std::string& file_name,
                        const std::string& compressor,
                        const std::vector<std::string>& lines) {
  if (compressor.empty()) {
    std::ofstream file(file_name);

    if (file.fail()) {
      return absl::InternalError("failed to open file");
    }

    for (const std::string& line : lines) {
      file << line << '\n';
    }

    return absl::OkStatus();
  }

  std::FILE* f =
      popen(absl::StrFormat("%s > %s", compressor, file_name).c_str(), "w");

  if (f == NULL) {
    return absl::InternalError("failed to open a sub-process");
  }

  // Writes data to the process.
  for (const std::string& line : lines) {
    if (std::fputs(line.c_str(), f) == EOF) {
      return absl::InternalError("failed to write to the process");
    }

    if (std::fputc('\n', f) == EOF) {
      return absl::InternalError("failed to write to the process");
    }
  }

  // Waits for the process to exit.
  {
    const int exit_status = pclose(f);

    if (exit_status != 0) {
      return absl::InternalError(absl::StrFormat(
          "process failed with non-zero exit code: %d", exit_status));
    }
  }

  return absl::OkStatus();
}

#endif
