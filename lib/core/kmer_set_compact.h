#ifndef CORE_KMER_SET_COMPACT_H_
#define CORE_KMER_SET_COMPACT_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <string>
#include <tuple>
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

// KmerSetCompact can store a set of kmers efficiently by using the binary form
// of SPSS. It is not possible to add or remove kmers from the structure, but it
// is possible to dump data to a file or load data from a file.
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

    return KmerSetCompact(spss, n_workers);
  }

  // Constructs a KmerSet.
  // If the structure was created with "canonical == true", "canonical" should
  // be true here as well.
  KmerSet<K, N, KeyType> ToKmerSet(bool canonical, int n_workers) const {
    return GetKmerSetFromSPSS<K, N, KeyType>(ToStrings(n_workers), canonical,
                                             n_workers);
  }

  // Dumps data to a file.
  // Dumped data can be read by Load().
  // Example:
  //   ... = Dump("foo.kmersetcompact.bz2", "bzip2", n_workers);
  //   ... = Dump("foo.kmersetcompact", "", n_workers);
  absl::Status Dump(const std::string& file_name, const std::string& compressor,
                    int n_workers) const {
    return WriteLines(file_name, compressor, ToStrings(n_workers));
  }

  // Loads data from a file.
  // Example:
  //   ... = Load("foo.kmersetcompact.bz2", "bzip2 -d", n_workers);
  //   ... = Load("foo.kmersetcompact", "", n_workers);
  static absl::StatusOr<KmerSetCompact> Load(const std::string& file_name,
                                             const std::string& decompressor,
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

    return KmerSetCompact(lines, n_workers);
  }

  // Returns the number of stored kmers.
  std::int64_t Size() const {
    std::int64_t size = 0;

    for (const std::vector<bool>& v : spss_binary_) {
      size += v.size() / 2 - K + 1;
    }

    return size;
  }

  // Constructs a structure like KmerSet, with selected buckets.
  // Each bucket is a vector of KeyType instead of a hash set, and they are
  // sorted.
  std::vector<std::vector<KeyType>> GetSampledKmerSet(
      const std::vector<int>& bucket_ids, bool canonical, int n_workers) const {
    std::vector<std::string> spss = ToStrings(n_workers);

    int n_buckets = bucket_ids.size();

    // bucket_ids[i] = j <-> map[j] = i.
    absl::flat_hash_map<int, int> map;
    for (int i = 0; i < n_buckets; i++) {
      map[bucket_ids[i]] = i;
    }

    std::vector<std::vector<KeyType>> buckets(n_buckets);

    {
      boost::asio::thread_pool pool(n_workers);
      std::vector<std::mutex> mus(n_buckets);

      for (const Range& range :
           Range(0, spss.size()).Split(n_workers * n_workers)) {
        boost::asio::post(pool, [&, range] {
          std::vector<std::vector<KeyType>> buf(n_buckets);

          for (std::int64_t i : range) {
            const std::string& s = spss[i];

            for (int j = 0; j < static_cast<int>(s.length()) - K + 1; j++) {
              Kmer<K> kmer(s.substr(j, K));

              if (canonical) kmer = kmer.Canonical();

              int bucket;
              KeyType key;
              std::tie(bucket, key) =
                  GetBucketAndKeyFromKmer<K, N, KeyType>(kmer);

              auto it = map.find(bucket);
              if (it == map.end()) continue;

              buf[it->second].push_back(key);
            }
          }

          // Moves from buffer.

          if (n_workers == 1) {
            buckets = std::move(buf);
          } else {
            int done_count = 0;
            std::vector<bool> done(n_buckets);

            while (done_count < n_buckets) {
              for (int i = 0; i < n_buckets; i++) {
                if (!done[i] && mus[i].try_lock()) {
                  buckets[i].insert(buckets[i].end(), buf[i].begin(),
                                    buf[i].end());

                  mus[i].unlock();
                  done[i] = true;
                  done_count += 1;
                }
              }
            }
          }
        });
      }

      pool.join();
    }

    // Sorts keys in each bucket.
    {
      boost::asio::thread_pool pool(n_workers);

      for (int i = 0; i < n_buckets; i++) {
        boost::asio::post(
            pool, [&, i] { std::sort(buckets[i].begin(), buckets[i].end()); });
      }

      pool.join();
    }

    return buckets;
  }

 private:
  explicit KmerSetCompact(std::vector<std::vector<bool>> spss_binary)
      : spss_binary_(std::move(spss_binary)){};

  KmerSetCompact(std::vector<std::string>& spss, int n_workers) {
    std::vector<std::vector<bool>> spss_binary(spss.size());

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, spss.size()).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (std::int64_t i : range) {
          const std::string& s = spss[i];

          spss_binary[i].resize(s.length() * 2);

          for (int j = 0; j < static_cast<int>(s.length()); j++) {
            switch (s[j]) {
              case 'A':
                break;
              case 'C':
                spss_binary[i][j * 2 + 1] = true;
                break;
              case 'G':
                spss_binary[i][j * 2] = true;
                break;
              case 'T':
                spss_binary[i][j * 2] = true;
                spss_binary[i][j * 2 + 1] = true;
                break;
              default:
                assert(false);
            }
          }
        }
      });
    }

    pool.join();

    spss_binary_ = std::move(spss_binary);
  }

  // Returns the SPSS as a vector of strings.
  std::vector<std::string> ToStrings(int n_workers) const {
    std::vector<std::string> spss(spss_binary_.size());

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, spss_binary_.size()).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (std::int64_t i : range) {
          spss[i].resize(spss_binary_[i].size() / 2);

          for (int j = 0; j < static_cast<int>(spss_binary_[i].size()) / 2;
               j++) {
            bool first = spss_binary_[i][j * 2];
            bool second = spss_binary_[i][j * 2 + 1];

            if (first) {
              if (second) {
                spss[i][j] = 'T';
              } else {
                spss[i][j] = 'G';
              }
            } else {
              if (second) {
                spss[i][j] = 'C';
              } else {
                spss[i][j] = 'A';
              }
            }
          }
        }
      });
    }

    pool.join();

    return spss;
  }

  std::vector<std::vector<bool>> spss_binary_;
};

#endif