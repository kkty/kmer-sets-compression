#ifndef KMER_SET_H_
#define KMER_SET_H_

#include <array>
#include <mutex>
#include <string>
#include <thread>
#include <random>
#include <vector>

#include "absl/strings/str_split.h"
#include "kmer.h"
#include "range.h"

// Use (1 << B) buckets internally.
template <int K, int B>
class KmerSet {
 public:
  // Find k-mers in a FASTQ file.
  void FromFASTQ(std::istream& is) {
    std::vector<std::string> lines;

    {
      std::array<std::string, 4> buf;

      // Reads 4 lines from "is" and pushes to "buf".
      auto consume = [&is, &buf]() {
        for (int i = 0; i < 4; i++) {
          if (!std::getline(is, buf[i])) return false;
        }

        return true;
      };

      while (consume() && buf[0].length() && buf[0][0] == '@') {
        lines.push_back(std::move(buf[1]));
      }
    }

    std::array<std::mutex, 1 << B> buckets_mutexes;

    std::vector<std::thread> threads;

    const auto ranges =
        SplitRange(0, lines.size(), std::thread::hardware_concurrency());

    for (const auto& range : ranges) {
      threads.emplace_back([&, range] {
        std::array<Bucket, 1 << B> buf;

        for (int64_t i = range.first; i < range.second; i++) {
          const auto& line = lines[i];
          std::vector<std::string> fragments = absl::StrSplit(line, "N");
          for (const auto& fragment : fragments) {
            for (int i = 0; i + K <= (int)fragment.length(); i++) {
              const Kmer<K> kmer(fragment.substr(i, K));

              int bucket;
              std::bitset<K * 2 - B> key;
              std::tie(bucket, key) = GetBucketAndKey(kmer);

              buf[bucket].insert(key);
            }
          }
        }

        std::vector<int> v;
        v.reserve(1 << B);
        for (int i = 0; i < (1 << B); i++) v.push_back(i);

        {
          std::default_random_engine e(range.first);
          std::shuffle(v.begin(), v.end(), e);
        }

        for (int i : v) {
          std::lock_guard _{buckets_mutexes[i]};
          auto& bucket = buckets_[i];
          bucket.reserve(bucket.size() + buf[i].size());
          for (const auto& key : buf[i]) {
            bucket.insert(key);
          }
        }
      });
    }

    for (std::thread& thread : threads) {
      thread.join();
    }
  }

  // Loads from dumped data.
  // Returns "true" if successful.
  bool Load(std::istream& is) {
    int k, b;
    is >> k >> b;

    if (k != K || b != B) {
      return false;
    }

    for (int i = 0; i < (1 << B); i++) {
      int64_t size;
      is >> size;
      buckets_[i].reserve(size);

      for (int64_t j = 0; j < size; j++) {
        std::string kmer_s;
        is >> kmer_s;
        std::bitset<K * 2 - B> kmer_b(kmer_s);
        buckets_[i].insert(std::move(kmer_b));
      }
    }

    return true;
  }

  void Dump(std::ostream& os) const {
    os << K << ' ' << B << '\n';

    for (int i = 0; i < (1 << B); i++) {
      os << buckets_[i].size() << '\n';
      for (const auto& key : buckets_[i]) {
        os << key.to_string() << '\n';
      }
    }
  }

  int Size() const {
    int sum = 0;
    for (const auto& bucket : buckets_) {
      sum += bucket.size();
    }
    return sum;
  }

  bool Contains(const Kmer<K>& kmer) const {
    int bucket;
    std::bitset<K * 2 - B> key;
    std::tie(bucket, key) = GetBucketAndKey(kmer);
    return buckets_[bucket].find(key) != buckets_[bucket].end();
  }

  Kmer<K> Sample() const {
    const int bucket = rand() % (1 << B);
    const int64_t n = rand() % buckets_[bucket].size();

    auto it = buckets_[bucket].begin();
    for (int64_t i = 0; i < n; i++) it++;

    const auto key = *it;

    std::bitset<K * 2> kmer;

    for (int i = 0; i < K * 2 - B; i++) kmer.set(i, key[i]);
    {
      std::bitset<B> bucket_b(bucket);
      for (int i = 0; i < B; i++) kmer.set(i + K * 2 - B, bucket_b[i]);
    }

    return Kmer<K>{kmer};
  }

 private:
  using Bucket = absl::flat_hash_set<std::bitset<K * 2 - B>>;
  std::array<Bucket, 1 << B> buckets_;

  static std::pair<int, std::bitset<K * 2 - B>> GetBucketAndKey(
      const Kmer<K>& kmer) {
    const auto& bits = kmer.Bits();
    const int key_bits = K * 2 - B;
    const int bucket = bits.to_ullong() >> key_bits;
    const std::bitset<K * 2 - B> key{bits.to_ullong() % (1ull << key_bits)};
    return {bucket, key};
  }
};

#endif
