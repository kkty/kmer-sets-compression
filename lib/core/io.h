#ifndef IO_H_
#define IO_H_

#include <fstream>
#include <istream>
#include <string>
#include <vector>

#include "boost/process.hpp"

std::vector<std::string> ReadLines(std::istream& is) {
  std::vector<std::string> lines;

  std::string s;
  while (std::getline(is, s)) {
    lines.push_back(std::move(s));
  }

  return lines;
}

// Opens a file and reads lines in it.
std::vector<std::string> ReadLines(const std::string& file_name) {
  std::ifstream file(file_name);
  return ReadLines(file);
}

// Opens a file, decompresses it, and reads lines in it.
// Example: ReadLines("foo.gz", "gzcat")
std::vector<std::string> ReadLines(const std::string& file_name,
                                   const std::string& decompressor) {
  boost::process::ipstream out;

  boost::process::child child(
      decompressor,
      boost::process::std_in<file_name, boost::process::std_out> out,
      boost::process::std_err > boost::process::null);

  std::vector<std::string> lines = ReadLines(out);

  child.wait();

  return lines;
}

void WriteLines(std::ostream& os, const std::vector<std::string>& lines) {
  for (const auto& line : lines) os << line << '\n';
}

void WriteLines(const std::string& file_name,
                const std::vector<std::string>& lines) {
  std::ofstream file(file_name);
  WriteLines(file, lines);
}

void WriteLines(const std::string& file_name, const std::string& compressor,
                const std::vector<std::string>& lines) {
  boost::process::opstream in;

  boost::process::child child(
      compressor, boost::process::std_in<in, boost::process::std_out> file_name,
      boost::process::std_err > boost::process::null);

  WriteLines(in, lines);

  in.close();
  in.pipe().close();

  child.wait();
}

#endif
