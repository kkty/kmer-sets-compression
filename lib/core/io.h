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

  os << std::flush;
}

}  // namespace internal

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

    return internal::ReadLines(file);
  }

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

    internal::WriteLines(file, lines);

    return absl::OkStatus();
  }

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
