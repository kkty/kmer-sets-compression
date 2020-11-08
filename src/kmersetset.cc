#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/status/status.h"
#include "boost/iostreams/filter/bzip2.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "boost/process.hpp"
#include "graph.h"
#include "kmer.h"
#include "kmer_counter.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"
#include "kmer_set_set.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for FASTQ files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");

int main(int argc, char** argv) {
  spdlog::set_default_logger(spdlog::stderr_color_mt("default"));

  // List of FASTQ file names.
  const std::vector<std::string> files = [&] {
    std::vector<std::string> v;

    for (const auto arg : absl::ParseCommandLine(argc, argv)) {
      v.push_back(arg);
    }

    // Program name.
    v.erase(v.begin());
    return v;
  }();

  const int K = 21;
  const int B = 6;

  const int n_datasets = files.size();

  std::vector<KmerSet<K, B>> kmer_sets(n_datasets);

  for (int i = 0; i < n_datasets; i++) {
    const std::string& file = files[i];

    spdlog::info("constructing kmer_counter for {}", file);

    KmerCounter<K, B> kmer_counter;

    {
      const std::string decompressor = absl::GetFlag(FLAGS_decompressor);

      if (decompressor != "") {
        boost::process::ipstream ipstream;
        boost::process::child child(decompressor + ' ' + file,
                                    boost::process::std_out > ipstream);

        const absl::Status status =
            kmer_counter.FromFASTQ(ipstream, absl::GetFlag(FLAGS_canonical));

        child.wait();

        if (!status.ok()) {
          spdlog::error("failed to parse FASTQ file");
          std::exit(1);
        }
      } else {
        std::ifstream is{file};

        const absl::Status status =
            kmer_counter.FromFASTQ(is, absl::GetFlag(FLAGS_canonical));

        if (!status.ok()) {
          spdlog::error("failed to parse FASTQ file");
          std::exit(1);
        }
      }
    }

    spdlog::info("constructed kmer_counter for {}", file);
    spdlog::info("constructing kmer_set for {}", file);

    const auto [kmer_set, cutoff_count] =
        kmer_counter.Set(absl::GetFlag(FLAGS_cutoff));

    spdlog::info("constructed kmer_set for {}", file);
    spdlog::info("cutoff_count = {}", cutoff_count);

    kmer_sets[i] = std::move(kmer_set);
  }

  int64_t kmer_sets_size = 0;

  for (int i = 0; i < n_datasets; i++) {
    int64_t kmer_set_size = kmer_sets[i].Size();
    spdlog::info("i = {}, kmer_set_size = {}", i, kmer_set_size);
    kmer_sets_size += kmer_set_size;
  }

  spdlog::info("kmer_sets_size = {}", kmer_sets_size);

  spdlog::info("constructing kmer_set_set");
  KmerSetSet<K, B> kmer_set_set(kmer_sets);
  spdlog::info("constructed kmer_set_set");

  spdlog::info("kmer_set_set.Size() = {}", kmer_set_set.Size());

  // Make sure that we can re-construct original k-mer sets.

  for (int i = 0; i < n_datasets; i++) {
    spdlog::info("kmer_sets[{}] == kmer_set_set.Get({}) = {}", i, i,
                 kmer_sets[i] == kmer_set_set.Get(i));
  }

  return 0;
}
