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

ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for FASTQ files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");
ABSL_FLAG(bool, fastcost, false, "use fast cost calculation");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(int, recursion, 1, "recursion limit for KmerSetSet");

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

  if (absl::GetFlag(FLAGS_debug)) spdlog::set_level(spdlog::level::debug);

  const int n_workers = absl::GetFlag(FLAGS_workers);

  const int K = 15;
  using KeyType = uint16_t;

  const int n_datasets = files.size();

  std::vector<KmerSet<K, KeyType>> kmer_sets(n_datasets);

  for (int i = 0; i < n_datasets; i++) {
    const std::string& file = files[i];

    spdlog::info("constructing kmer_counter for {}", file);

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

    spdlog::info("constructed kmer_counter for {}", file);
    spdlog::info("constructing kmer_set for {}", file);

    const auto [kmer_set, cutoff_count] =
        kmer_counter.Set(absl::GetFlag(FLAGS_cutoff), n_workers);

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

  const auto cost_function = [fastcost = absl::GetFlag(FLAGS_fastcost),
                              canonical = absl::GetFlag(FLAGS_canonical)](
                                 const KmerSet<K, KeyType>& lhs,
                                 const KmerSet<K, KeyType>& rhs,
                                 int n_workers) {
    int64_t cost = 0;

    if (fastcost) {
      cost += Sub(lhs, rhs, n_workers).Size();
      cost += Sub(rhs, lhs, n_workers).Size();
    } else if (canonical) {
      for (const std::string& unitig :
           GetUnitigsCanonical(Sub(lhs, rhs, n_workers), n_workers))
        cost += unitig.length();
      for (const std::string& unitig :
           GetUnitigsCanonical(Sub(rhs, lhs, n_workers), n_workers))
        cost += unitig.length();
    } else {
      for (const std::string& unitig :
           GetUnitigs(Sub(lhs, rhs, n_workers), n_workers))
        cost += unitig.length();
      for (const std::string& unitig :
           GetUnitigs(Sub(rhs, lhs, n_workers), n_workers))
        cost += unitig.length();
    }

    return cost;
  };

  {
    int64_t total_cost = 0;
    for (int i = 0; i < n_datasets; i++) {
      int64_t cost =
          cost_function(KmerSet<K, KeyType>(), kmer_sets[i], n_workers);
      spdlog::info("i = {}, cost = {}", i, cost);
      total_cost += cost;
    }
    spdlog::info("total_cost = {}", total_cost);
  }

  spdlog::info("constructing kmer_set_set");

  KmerSetSet<K, KeyType, decltype(cost_function)> kmer_set_set(
      std::move(kmer_sets), absl::GetFlag(FLAGS_recursion), cost_function,
      n_workers);
  spdlog::info("constructed kmer_set_set");

  spdlog::info("kmer_set_set.Size() = {}", kmer_set_set.Size());
  spdlog::info("kmer_set_set.Depth() = {}", kmer_set_set.Depth());
  spdlog::info("kmer_set_set.Cost() = {}", kmer_set_set.Cost());

  // Make sure that we can re-construct original k-mer sets.

  for (int i = 0; i < n_datasets; i++) {
    spdlog::info(
        "Equals(kmer_sets[{}], kmer_set_set.Get({})) = {}", i, i,
        Equals(kmer_sets[i], kmer_set_set.Get(i, n_workers), n_workers));
  }

  return 0;
}
