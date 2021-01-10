#ifndef CORE_KMER_SET_COMPACT_H_
#define CORE_KMER_SET_COMPACT_H_

#include <algorithm>
#include <atomic>
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
#include "streamvbyte.h"

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

    return KmerSetCompact(spss);
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

    return KmerSetCompact(lines);
  }

  // Returns the number of stored kmers.
  std::int64_t Size(int n_workers) const {
    std::vector<std::uint32_t> lengths = GetLengths(n_workers);

    std::atomic_int64_t size = 0;

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range : Range(0, n_).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        std::int64_t buf = 0;

        for (std::int64_t i : range) {
          buf += lengths[i] - K + 1;
        }

        size += buf;
      });
    }

    pool.join();

    return size;
  }

  // Returns the sum of each string's length.
  std::int64_t Weight(std::int64_t n_workers) const {
    std::atomic_int64_t weight = 0;

    std::vector<std::uint32_t> lengths = GetLengths(n_workers);

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range :
         Range(0, lengths.size()).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        std::int64_t buf = 0;

        for (std::int64_t i : range) {
          buf += lengths[i];
        }

        weight += buf;
      });
    }

    pool.join();

    return static_cast<std::int64_t>(weight);
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
  KmerSetCompact(std::vector<std::string>& spss) {
    n_ = spss.size();

    // lengths[i] stored the length of spss[i], subtracted by K.
    std::vector<std::uint32_t> lengths(n_);

    // positions[i] stores the starting position in data_ for the ith string.
    std::vector<std::int64_t> positions(n_);

    // Calculates positions, lengths, and size of data_.
    {
      std::int64_t size = 0;

      for (std::int64_t i = 0; i < n_; i++) {
        positions[i] = size;

        std::uint32_t length = spss[i].length();
        lengths[i] = length - K;
        size += length * 2;
      }

      data_.resize(size);
    }

    for (std::int64_t i = 0; i < n_; i++) {
      const std::string& s = spss[i];

      std::int64_t position = positions[i];
      int length = lengths[i] + K;

      for (int j = 0; j < length; j++) {
        switch (s[j]) {
          case 'A':
            break;
          case 'C':
            data_[position + j * 2 + 1] = true;
            break;
          case 'G':
            data_[position + j * 2] = true;
            break;
          case 'T':
            data_[position + j * 2] = true;
            data_[position + j * 2 + 1] = true;

            break;
          default:
            assert(false);
        }
      }
    }

    {
      lengths_compressed_.resize(streamvbyte_max_compressedbytes(n_));

      std::int64_t size = streamvbyte_encode_0124(lengths.data(), n_,
                                                  lengths_compressed_.data());

      lengths_compressed_.resize(size);
      lengths_compressed_.shrink_to_fit();
    }
  }

  std::vector<std::uint32_t> GetLengths(int n_workers) const {
    std::vector<std::uint32_t> lengths(n_);

    streamvbyte_decode_0124(lengths_compressed_.data(), lengths.data(), n_);

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range : Range(0, n_).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (std::int64_t i : range) {
          lengths[i] += K;
        }
      });
    }

    pool.join();

    return lengths;
  }

  // Returns the SPSS as a vector of strings.
  std::vector<std::string> ToStrings(int n_workers) const {
    std::vector<std::uint32_t> lengths = GetLengths(n_workers);

    std::vector<std::int64_t> positions(n_);

    for (std::int64_t i = 0; i < n_ - 1; i++) {
      positions[i + 1] = positions[i] + lengths[i] * 2;
    }

    std::vector<std::string> spss(n_);

    boost::asio::thread_pool pool(n_workers);

    for (const Range& range : Range(0, n_).Split(n_workers * n_workers)) {
      boost::asio::post(pool, [&, range] {
        for (std::int64_t i : range) {
          const int length = lengths[i];
          const std::int64_t position = positions[i];

          spss[i].resize(length);

          for (int j = 0; j < length; j++) {
            const bool first = data_[position + 2 * j];
            const bool second = data_[position + 2 * j + 1];

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

  // Number of stored strings.
  std::int64_t n_;

  // Lengths of each string.
  // Its size should be equal to n_.
  std::vector<std::uint8_t> lengths_compressed_;

  // All the strings are saved to this one vector. It is more memory efficient
  // than std::vector<std::vector<bool>>.
  std::vector<bool> data_;
};

#endif