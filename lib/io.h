#ifndef IO_H_
#define IO_H_

#include <string>
#include <tuple>
#include <thread>
#include <vector>

#include "absl/status/statusor.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compressed.h"
#include "spdlog/spdlog.h"

template <int K, typename KeyType>
KmerSet<K, KeyType> GetKmerSetFromCompressedKmersFile(
    const std::string& file_name, const std::string& decompressor,
    bool canonical, int n_workers) {
  spdlog::info("constructing kmer_set_compressed");

  const KmerSetCompressed<K, KeyType> kmer_set_compressed =
      KmerSetCompressed<K, KeyType>::Load(file_name, decompressor);

  spdlog::info("constructed kmer_set_compressed");

  spdlog::info("constructing kmer_set");

  const KmerSet<K, KeyType> kmer_set =
      kmer_set_compressed.ToKmerSet(canonical, n_workers);

  spdlog::info("constructed kmer_set");

  return kmer_set;
}

template <int K, typename KeyType>
KmerSet<K, KeyType> GetKmerSetFromFASTQFile(const std::string& file_name,
                                            const std::string& decompressor,
                                            bool canonical, int cutoff,
                                            int n_workers) {
  spdlog::info("constructing kmer_counter");

  const KmerCounter<K, KeyType> kmer_counter = [&] {
    absl::StatusOr<KmerCounter<K, KeyType>> statusor =
        decompressor != "" ? KmerCounter<K, KeyType>::FromFASTQ(
                                 file_name, decompressor, canonical, n_workers)
                           : KmerCounter<K, KeyType>::FromFASTQ(
                                 file_name, canonical, n_workers);

    if (!statusor.ok()) {
      spdlog::error("failed to parse FASTQ file: {}",
                    statusor.status().ToString());
      std::exit(1);
    }

    return *statusor;
  }();

  spdlog::info("constructed kmer_counter");

  spdlog::info("constructing kmer_set");

  KmerSet<K, KeyType> kmer_set;
  int64_t cutoff_count;
  std::tie(kmer_set, cutoff_count) = kmer_counter.ToKmerSet(cutoff, n_workers);

  spdlog::info("constructed kmer_set");

  return kmer_set;
}

// Returns a list of reads in a FASTA file.
std::vector<std::string> ReadFASTAFile(const std::string& file_name,
                                       const std::string& decompressor,
                                       int n_workers) {
  std::vector<std::string> lines = ReadLines(file_name, decompressor);

  const int n_lines = lines.size();

  std::vector<std::string> reads(n_lines / 2);

  {
    std::vector<std::thread> threads;

    for (const Range& range : Range(0, n_lines).Split(n_workers)) {
      threads.emplace_back([&, range] {
        range.ForEach([&](int64_t i) {
          const std::string& line = lines[i];
          if (i % 2 == 0) {
            if (line[0] != '>') {
              spdlog::error("invalid file");
              std::exit(1);
            }
            std::string().swap(lines[i]);
          } else {
            reads[i / 2] = std::move(lines[i]);
          }
        });
      });
    }

    for (std::thread& t : threads) t.join();
  }

  return reads;
}

#endif