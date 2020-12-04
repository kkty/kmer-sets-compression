#ifndef CORE_IO_H_
#define CORE_IO_H_

#include <fstream>
#include <istream>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "boost/process.hpp"

namespace internal {

std::vector<std::string> ReadLines(std::istream& is) {
  std::vector<std::string> lines;

  std::string s;
  while (std::getline(is, s)) {
    lines.push_back(s);
  }

  return lines;
}

void WriteLines(std::ostream& os, const std::vector<std::string>& lines) {
  for (const std::string& line : lines) {
    os << line << '\n';
  }
}

}  // namespace internal

// Opens a file and reads lines in it.
absl::StatusOr<std::vector<std::string>> ReadLines(
    const std::string& file_name) {
  std::ifstream file(file_name);

  if (file.fail()) {
    return absl::InternalError("failed to open file");
  }

  return internal::ReadLines(file);
}

// Opens a file, decompresses it, and reads lines in it.
// Example: ReadLines("foo.gz", "gzcat")
absl::StatusOr<std::vector<std::string>> ReadLines(
    const std::string& file_name, const std::string& decompressor) {
  boost::process::ipstream out;

  boost::process::child child(
      decompressor,
      boost::process::std_in<file_name, boost::process::std_out> out,
      boost::process::std_err > boost::process::null);

  std::vector<std::string> lines = internal::ReadLines(out);

  child.wait();

  const int exit_code = child.exit_code();

  if (exit_code != 0) {
    return absl::InternalError(absl::StrFormat(
        "decompressor failed with non-zero exit code: %d", exit_code));
  }

  return lines;
}

absl::Status WriteLines(const std::string& file_name,
                        const std::vector<std::string>& lines) {
  std::ofstream file(file_name);

  if (file.fail()) {
    return absl::InternalError("failed to open file");
  }

  internal::WriteLines(file, lines);

  return absl::OkStatus();
}

absl::Status WriteLines(const std::string& file_name,
                        const std::string& compressor,
                        const std::vector<std::string>& lines) {
  boost::process::opstream in;

  boost::process::child child(
      compressor, boost::process::std_in<in, boost::process::std_out> file_name,
      boost::process::std_err > boost::process::null);

  internal::WriteLines(in, lines);

  in.close();
  in.pipe().close();

  child.wait();

  const int exit_code = child.exit_code();

  if (exit_code != 0) {
    return absl::InternalError(absl::StrFormat(
        "compressor failed with non-zero exit code: %d", exit_code));
  }

  return absl::OkStatus();
}

#endif
