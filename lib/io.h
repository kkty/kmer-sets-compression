#ifndef IO_H_
#define IO_H_

#include <atomic>
#include <cstdint>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "core/kmer_counter.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "spdlog/spdlog.h"

template <int K, int N, typename KeyType>
absl::StatusOr<KmerSet<K, N, KeyType>> GetKmerSetFromCompressedKmersFile(
    const std::string& file_name, const std::string& decompressor,
    bool canonical, int n_workers) {
  KmerSetCompact<K, N, KeyType> kmer_set_compact;

  {
    spdlog::info("constructing kmer_set_compact");

    absl::StatusOr<KmerSetCompact<K, N, KeyType>> statusor =
        KmerSetCompact<K, N, KeyType>::Load(file_name, decompressor);

    if (!statusor.ok()) {
      return statusor.status();
    }

    kmer_set_compact = std::move(statusor).value();

    spdlog::info("constructed kmer_set_compact");
  }

  spdlog::info("constructing kmer_set");

  const KmerSet<K, N, KeyType> kmer_set =
      kmer_set_compact.ToKmerSet(canonical, n_workers);

  spdlog::info("constructed kmer_set");

  return kmer_set;
}

// Returns a list of reads in a FASTA file.
absl::StatusOr<std::vector<std::string>> ReadFASTAFile(
    const std::string& file_name, const std::string& decompressor,
    int n_workers) {
  std::vector<std::string> lines;

  {
    absl::StatusOr<std::vector<std::string>> statusor =
        ReadLines(file_name, decompressor);

    if (!statusor.ok()) {
      return statusor.status();
    }

    lines = std::move(statusor).value();
  }

  const int n_lines = lines.size();

  if (n_lines % 2 != 0) {
    return absl::FailedPreconditionError(
        "FASTA files should contain an even number of lines");
  }

  std::vector<std::string> reads(n_lines / 2);

  {
    std::vector<std::thread> threads;
    std::atomic_bool is_valid = true;

    for (const Range& range : Range(0, n_lines).Split(n_workers)) {
      threads.emplace_back([&, range] {
        for (std::int64_t i : range) {
          const std::string& line = lines[i];
          if (i % 2 == 0) {
            if (line[0] != '>') {
              is_valid = false;
            }

            std::string().swap(lines[i]);
          } else {
            reads[i / 2] = std::move(lines[i]);
          }
        }
      });
    }

    for (std::thread& t : threads) t.join();

    if (!is_valid) {
      return absl::FailedPreconditionError("invalid FASTA file");
    }
  }

  return reads;
}

#endif