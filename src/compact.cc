#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/status/status.h"
#include "boost/iostreams/filter/bzip2.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/process.hpp"
#include "kmer.h"
#include "kmer_counter.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");
ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for FASTQ files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, workers, 1, "number of workers");

int main(int argc, char** argv) {
  std::vector<std::string> files = [&] {
    std::vector<std::string> v;
    for (const char* arg : absl::ParseCommandLine(argc, argv)) {
      v.push_back(arg);
    }
    v.erase(v.begin());
    return v;
  }();

  if (absl::GetFlag(FLAGS_debug)) spdlog::set_level(spdlog::level::debug);

  std::srand(std::time(nullptr));
  std::ios_base::sync_with_stdio(false);

  const int K = 21;
  const int n_workers = absl::GetFlag(FLAGS_workers);
  using KeyType = uint32_t;

  for (const std::string& file : files) {
    spdlog::info("file = {}", file);

    spdlog::info("constructing kmer_counter");

    KmerCounter<K, KeyType> kmer_counter;

    {
      const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

      if (decompressor != "") {
        boost::process::ipstream ipstream;
        boost::process::child child(decompressor + ' ' + file,
                                    boost::process::std_out > ipstream);

        const absl::Status status = kmer_counter.FromFASTQ(
            ipstream, absl::GetFlag(FLAGS_canonical), n_workers);

        child.wait();

        if (!status.ok()) {
          spdlog::error("failed to parse FASTQ file");
          std::exit(1);
        }
      } else {
        std::ifstream is{file};

        const absl::Status status = kmer_counter.FromFASTQ(
            is, absl::GetFlag(FLAGS_canonical), n_workers);

        if (!status.ok()) {
          spdlog::error("failed to parse FASTQ file");
          std::exit(1);
        }
      }
    }

    spdlog::info("constructed kmer_counter");
    spdlog::info("kmer_counter.Size() = {}", kmer_counter.Size());

    spdlog::info("constructing kmer_set");

    const auto [kmer_set, cutoff_count] =
        kmer_counter.Set(absl::GetFlag(FLAGS_cutoff), n_workers);

    spdlog::info("kmer_set.Size() = {}", kmer_set.Size());

    spdlog::info("constructed kmer_set");
    spdlog::info("cutoff_count = {}", cutoff_count);

    const KmerSetCompact kmer_set_compact(
        kmer_set, absl::GetFlag(FLAGS_canonical), n_workers);

    spdlog::info("kmer_set_compact.Size() = {}", kmer_set_compact.Size());

    spdlog::info("kmer_set_compact.Find().size() = {}",
                 kmer_set_compact.Find(n_workers).size());

    {
      std::vector<Kmer<K>> original = kmer_set.Find(n_workers);
      std::vector<Kmer<K>> from_compact = kmer_set_compact.Find(n_workers);

      std::sort(original.begin(), original.end());
      std::sort(from_compact.begin(), from_compact.end());

      spdlog::info("original == from_compact = {}", original == from_compact);
    }
  }
}
