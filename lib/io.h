#ifndef IO_H_
#define IO_H_

#include <cstdint>
#include <filesystem>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <utility>

#include "absl/random/random.h"
#include "absl/status/statusor.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "spdlog/spdlog.h"

// Loads a KmerSetCompact from a file and constructs a KmerSet from it.
template <int K, int N, typename KeyType>
absl::StatusOr<KmerSet<K, N, KeyType>> GetKmerSetFromFile(
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

// Creates a path to a temporary file on initialization and removes it on
// deconstruction.
class TemporaryFile {
 public:
  TemporaryFile() {
    std::stringstream ss;
    absl::InsecureBitGen bitgen;

    ss << std::hex
       << absl::Uniform<std::uint64_t>(
              absl::IntervalClosed, bitgen, 0,
              std::numeric_limits<std::uint64_t>::max());

    name_ = (std::filesystem::temp_directory_path() / ss.str()).string();
  }

  ~TemporaryFile() { std::filesystem::remove(name_); }

  // Returns the file path.
  std::string Name() const { return name_; }

 private:
  std::string name_;
};

// Creates a path to a temporary directory on initialization and removes it
// on deconstruction.
class TemporaryDirectory {
 public:
  TemporaryDirectory() {
    std::stringstream ss;
    absl::InsecureBitGen bitgen;

    ss << std::hex
       << absl::Uniform<std::uint64_t>(
              absl::IntervalClosed, bitgen, 0,
              std::numeric_limits<std::uint64_t>::max());

    name_ = (std::filesystem::temp_directory_path() / ss.str()).string();
  }

  ~TemporaryDirectory() { std::filesystem::remove_all(name_); }

  // Returns the directory path.
  std::string Name() const { return name_; }

 private:
  std::string name_;
};

#endif