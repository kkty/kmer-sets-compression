#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "kmer.h"
#include "kmer_counter.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");

int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);

  const std::vector<std::string> files = [&] {
    std::vector<std::string> v;
    for (const auto arg : absl::ParseCommandLine(argc, argv)) {
      v.push_back(arg);
    }
    v.erase(v.begin());
    return v;
  }();

  if (absl::GetFlag(FLAGS_debug)) spdlog::set_level(spdlog::level::debug);

  const int K = 21;
  const int B = 6;

  const int n_datasets = files.size();

  std::vector<KmerSet<K, B>> kmer_sets(n_datasets);

  for (int i = 0; i < n_datasets; i++) {
    const std::string& file = files[i];
    std::ifstream is{file};

    spdlog::info("constructing kmer_counter for {}", file);

    KmerCounter<K, B> kmer_counter;
    kmer_counter.FromFASTQ(is);

    spdlog::info("constructed kmer_counter for {}", file);
    spdlog::info("constructing kmer_set for {}", file);

    const auto [kmer_set, cutoff_count] =
        kmer_counter.Set(absl::GetFlag(FLAGS_cutoff));

    spdlog::info("constructed kmer_set for {}", file);

    kmer_sets[i] = std::move(kmer_set);
  }

  for (int i = 0; i < n_datasets; i++) {
    spdlog::info("kmer_sets[{}].Size() = {}", i, kmer_sets[i].Size());
  }
}
