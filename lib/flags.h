#ifndef FLAGS_H_
#define FLAGS_H_

#include <string>
#include <vector>

#include "absl/flags/parse.h"

// Parses the command line flags.
// The positional arguments are returned as a vector of strings, with the
// program name (argv[0]) stripped off.
std::vector<std::string> ParseFlags(int argc, char* argv[]) {
  std::vector<std::string> v;

  for (const char* arg : absl::ParseCommandLine(argc, argv)) {
    v.push_back(arg);
  }

  v.erase(v.begin());

  return v;
}

// Returns a flag message shared by multiple executables.
std::string GetFlagMessage(const std::string& s) {
  if (s == "k") {
    return "the length of k-mers";
  }

  if (s == "debug") {
    return "enable debugging messages";
  }

  if (s == "compressor") {
    return "an program to compress output files; e.g., \"bzip2\" for bzip2, "
           "\"\" for no compression";
  }

  if (s == "decompressor") {
    return "an program to decompress input files; e.g., \"bzip2 -d\" for "
           "bzip2, \"\" for no decompression";
  }

  if (s == "workers") {
    return "number of threads to use";
  }

  if (s == "canonical") {
    return "set this flag when handling canonical k-mers";
  }

  return "";
}

#endif