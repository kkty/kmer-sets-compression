#ifndef CORE_KMER_SET_COMPACT_H_
#define CORE_KMER_SET_COMPACT_H_

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "boost/asio/post.hpp"
#include "boost/asio/thread_pool.hpp"
#include "core/kmer.h"
#include "core/kmer_set.h"
#include "core/range.h"
#include "core/spss.h"

// KmerSetCompact can store a set of kmers efficiently by using SPSS.
// It is not possible to add or remove kmers from the structure, but it is
// possible to dump data to a file or load data from a file.
template <int K, int N, typename KeyType>
class KmerSetCompact {
 public:
  KmerSetCompact() = default;

  // Constructs KmerSetCompact from a kmer set. If the kmer set only contains
  // canonical kmers, "canonical" should be true for better compression.
  static KmerSetCompact FromKmerSet(const KmerSet<K, N, KeyType>& kmer_set,
                                    bool canonical, bool fast, int n_workers) {
    std::vector<std::string> spss;

    if (canonical) {
      spss = GetSPSSCanonical<K, N, KeyType>(kmer_set, fast, n_workers);
    } else {
      spss = GetSPSS<K, N, KeyType>(kmer_set, n_workers);
    }

    return KmerSetCompact(std::move(spss));
  }

  // Constructs a KmerSet.
  // If the structure was created with "canonical == true", "canonical" should
  // be true here as well.
  KmerSet<K, N, KeyType> ToKmerSet(bool canonical, int n_workers) const {
    return GetKmerSetFromSPSS<K, N, KeyType>(spss_, canonical, n_workers);
  }

  // Constructs a list of kmers.
  // If the structure was created with "canonical == true", "canonical" should
  // be true here as well.
  std::vector<Kmer<K>> ToKmers(bool canonical, int n_workers) const {
    return GetKmersFromSPSS<K>(spss_, canonical, n_workers);
  }

  // Dumps data to a file.
  // Dumped data can be read by Load().
  // Example:
  //   ... = Dump("foo.kmersetcompact.bz2", "bzip2");
  //   ... = Dump("foo.kmersetcompact", "");
  absl::Status Dump(const std::string& file_name,
                    const std::string& compressor) const {
    return WriteLines(file_name, compressor, spss_);
  }

  // Dumps data to a vector of strings.
  // Dumped data can be read by Load().
  // Each string does not contain whitespaces or line breaks.
  std::vector<std::string> Dump() const { return spss_; }

  // Loads data from a file.
  // Example:
  //   ... = Load("foo.kmersetcompact.bz2", "bzip2 -d");
  //   ... = Load("foo.kmersetcompact", "");
  static absl::StatusOr<KmerSetCompact> Load(const std::string& file_name,
                                             const std::string& decompressor) {
    std::vector<std::string> lines;

    {
      absl::StatusOr<std::vector<std::string>> statusor =
          ReadLines(file_name, decompressor);

      if (!statusor.ok()) {
        return statusor.status();
      }

      lines = std::move(statusor).value();
    }

    return KmerSetCompact(std::move(lines));
  }

  // Loads data from a vector of strings.
  static KmerSetCompact Load(std::vector<std::string> v) {
    return KmerSetCompact(std::move(v));
  }

  // Returns the total length of stored strings.
  // If the data is dumped to a file without compression, the estimated file
  // size is Size() bytes.
  std::int64_t Size(int n_workers) const {
    std::vector<std::thread> threads;
    std::atomic_int64_t size = 0;

    for (const Range& range : Range(0, spss_.size()).Split(n_workers)) {
      threads.emplace_back([&, range] {
        for (int i : range) {
          size += spss_[i].length();
        }
      });
    }

    for (std::thread& t : threads) t.join();

    return size;
  }

 private:
  explicit KmerSetCompact(std::vector<std::string> spss)
      : spss_(std::move(spss)){};

  std::vector<std::string> spss_;
};

#endif