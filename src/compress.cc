#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "boost/iostreams/filter/bzip2.hpp"
#include "boost/iostreams/filtering_stream.hpp"
#include "graph.h"
#include "kmer.h"
#include "kmer_counter.h"
#include "kmer_set.h"
#include "kmer_set_compact.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(bool, bzip2, false, "accept bzip2-ed FASTQ files");
ABSL_FLAG(bool, canonical, false, "count canonical k-mers");
ABSL_FLAG(int, cutoff, 1, "cut off threshold");

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
  const int B = 6;

  const int n_datasets = files.size();

  std::vector<KmerSet<K, B>> kmer_sets(n_datasets);

  for (int i = 0; i < n_datasets; i++) {
    const std::string& file = files[i];

    spdlog::info("constructing kmer_counter for {}", file);

    KmerCounter<K, B> kmer_counter;

    if (absl::GetFlag(FLAGS_bzip2)) {
      std::ifstream is(file, std::ios_base::binary);
      boost::iostreams::filtering_istream f_is;
      f_is.push(boost::iostreams::bzip2_decompressor());
      f_is.push(is);
      kmer_counter.FromFASTQ(f_is, absl::GetFlag(FLAGS_canonical));
    } else {
      std::ifstream is{file};
      kmer_counter.FromFASTQ(is, absl::GetFlag(FLAGS_canonical));
    }

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

  // kmer_sets[n_datasets] is an empty set.
  kmer_sets.push_back(KmerSet<K, B>{});

  // Returns the required size to represent the difference between s1 and s2.
  const auto get_diff = [&](const KmerSet<K, B>& s1, const KmerSet<K, B>& s2) {
    // Returns the size after compaction.
    const auto get_size = [&](KmerSet<K, B> s) {
      const std::vector<std::string> unitigs = GetUnitigs(s);
      int64_t size = 0;
      for (const std::string& unitig : unitigs) size += unitig.size();
      return size;
    };

    return get_size(s1 - s2) + get_size(s2 - s1);
  };

  for (int i = 0; i < n_datasets; i++) {
    const int64_t diff = get_diff(kmer_sets[n_datasets], kmer_sets[i]);
    spdlog::info("i = {}, diff = {}", i, diff);
  }

  {
    for (int i = 0; i < n_datasets; i++) {
      for (int j = i + 1; j < n_datasets; j++) {
        spdlog::debug("calculating intersection: i = {}, j = {}", i, j);
        KmerSet<K, B> intersection =
            kmer_sets[i] - (kmer_sets[i] - kmer_sets[j]);
        spdlog::debug("calculated intersection: i = {}, j = {}", i, j);
        kmer_sets.push_back(std::move(intersection));
      }
    }
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
  const auto [cost, tree] = g.MST(n_datasets);

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
