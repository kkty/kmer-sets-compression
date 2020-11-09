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
#include "graph.h"
#include "kmer.h"
#include "kmer_counter.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"
#include "kmer_set_set.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "",
          "specify decompressor for FASTQ files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");
ABSL_FLAG(bool, fastdiff, false, "use fast diff calculation");
ABSL_FLAG(int, workers, 1, "number of workers");

int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);

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

  const int K = 21;
  const int n_workers = absl::GetFlag(FLAGS_workers);
  using KeyType = uint32_t;

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

    kmer_sets[i] = std::move(kmer_set);
  }

  for (int i = 0; i < n_datasets; i++) {
    spdlog::info("kmer_sets[{}].Size() = {}", i, kmer_sets[i].Size());
  }

  // kmer_sets[n_datasets] is an empty set.
  kmer_sets.push_back(KmerSet<K, KeyType>{});

  // Returns the required size to represent the difference between s1 and s2.
  const auto get_diff = [&](const KmerSet<K, KeyType>& s1,
                            const KmerSet<K, KeyType>& s2) {
    if (absl::GetFlag(FLAGS_fastdiff)) {
      return Sub(s1, s2, n_workers).Size() + Sub(s2, s1, n_workers).Size();
    }

    // Returns the size after compaction.
    const auto get_size = [&](const KmerSet<K, KeyType>& s) {
      const std::vector<std::string> unitigs = absl::GetFlag(FLAGS_canonical)
                                                   ? GetUnitigsCanonical(s, n_workers)
                                                   : GetUnitigs(s, n_workers);
      int64_t size = 0;
      for (const std::string& unitig : unitigs) size += unitig.size();
      return size;
    };

    return get_size(Sub(s1, s2, n_workers)) + get_size(Sub(s2, s1, n_workers));
  };

  for (int i = 0; i < n_datasets; i++) {
    const int64_t diff = get_diff(kmer_sets[n_datasets], kmer_sets[i]);
    spdlog::info("i = {}, diff = {}", i, diff);
  }

  BidirectionalGraph g;
  for (int i = 0; i < (int)kmer_sets.size(); i++) {
    for (int j = i + 1; j < (int)kmer_sets.size(); j++) {
      const int64_t diff = get_diff(kmer_sets[i], kmer_sets[j]);
      spdlog::info("i = {}, j = {}, diff = {}", i, j, diff);
      g.AddEdge(i, j, get_diff(kmer_sets[i], kmer_sets[j]));
    }
  }

  // Constructs the MST whose root is an empty set.
  auto [cost, tree] = g.MST(n_datasets);

  spdlog::info("cost = {}", cost);

  // BFS from root.
  {
    std::queue<int64_t> queue;
    queue.push(tree.Root());
    while (!queue.empty()) {
      int64_t n = queue.front();
      queue.pop();
      for (int64_t child : tree.Children(n)) {
        queue.push(child);
        spdlog::info("tree.Parent({}) = {}", child, n);
      }
    }
  }
}
