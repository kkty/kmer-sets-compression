#ifndef IO_H_
#define IO_H_

#include <istream>
#include <string>
#include <vector>

#include "boost/process.hpp"

std::vector<std::string> ReadLines(std::istream& is) {
  std::vector<std::string> lines;

  std::string s;
  while (std::getline(is, s)) {
    lines.push_back(s);
  }

  return lines;
}

// Opens a file and reads lines in it.
std::vector<std::string> ReadLines(const std::string& file_name) {
  std::ifstream is(file_name);

  return ReadLines(is);
}

// Opens a file, decompresses it, and reads lines in it.
// Example: ReadLines("foo.gz", "gzcat")
std::vector<std::string> ReadLines(const std::string& file_name,
                                   const std::string& decompressor) {
  boost::process::ipstream ipstream;
  boost::process::child child(decompressor + ' ' + file_name,
                              boost::process::std_out > ipstream);

  std::vector<std::string> lines = ReadLines(ipstream);

  child.wait();

  return lines;
}

#endif
