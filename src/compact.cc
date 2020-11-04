#include <algorithm>
#include <cstdlib>
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
  absl::ParseCommandLine(argc, argv);

  if (absl::GetFlag(FLAGS_debug)) spdlog::set_level(spdlog::level::debug);

  std::srand(std::time(nullptr));
  std::ios_base::sync_with_stdio(false);

  const int K = 21;
  const int B = 6;

  spdlog::info("constructing kmer_counter");

  KmerCounter<K, B> kmer_counter;
  kmer_counter.FromFASTQ(std::cin);

  spdlog::info("constructed kmer_counter");
  spdlog::info("kmer_counter.Size() = {}", kmer_counter.Size());

  spdlog::info("constructing kmer_set");

  const auto [kmer_set, cutoff_count] =
      kmer_counter.Set(absl::GetFlag(FLAGS_cutoff));

  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());

  spdlog::info("constructed kmer_set");
  spdlog::info("cutoff_count = {}", cutoff_count);

  const KmerSetCompact kmer_set_compact{kmer_set};

  spdlog::info("kmer_set_compact.Size() = {}", kmer_set_compact.Size());

  spdlog::info("kmer_set_compact.Find().size() = {}",
               kmer_set_compact.Find().size());

  {
    std::vector<Kmer<K>> original = kmer_set.Find();
    std::vector<Kmer<K>> from_compact = kmer_set_compact.Find();

    std::sort(original.begin(), original.end());
    std::sort(from_compact.begin(), from_compact.end());

    spdlog::info("original == from_compact = {}", original == from_compact);
  }
}
