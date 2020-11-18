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

#endif