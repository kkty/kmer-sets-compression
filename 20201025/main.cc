#include <algorithm>
#include <array>
#include <atomic>
#include <bitset>
#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/strings/str_split.h"
#include "absl/time/time.h"
#include "omp.h"
#include "spdlog/spdlog.h"

// Use (1 << B) buckets internally.
template <int K = 31, int B = 6>
class KmerSet {
 public:
  KmerSet(std::istream& is, const int64_t lines_approximate,
          const int64_t kmers_approximate) {
    std::atomic_bool done = false;
    std::vector<std::thread> threads;
    std::queue<std::string> queue;
    std::mutex queue_mutex;
    std::array<std::mutex, 1 << B> buckets_mutexes;

    for (auto& bucket : buckets_) {
      bucket.reserve(kmers_approximate >> B);
    }

    std::vector<std::string> lines;
    lines.reserve(lines_approximate);

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

    const int64_t n_lines = lines.size();

#pragma omp parallel for
    for (int i = 0; i < n_lines; i++) {
      const auto& line = lines[i];
      std::vector<std::string> fragments = absl::StrSplit(line, "N");
      for (const auto& fragment : fragments) {
        for (int i = 0; i + K <= fragment.length(); i++) {
          const std::string kmer_s = fragment.substr(i, K);
          std::bitset<K * 2> kmer_b;

          for (int j = 0; j < K; j++) {
            switch (kmer_s[j]) {
              case 'A':
                break;
              case 'C':
                kmer_b.set(j * 2);
                break;
              case 'G':
                kmer_b.set(j * 2 + 1);
                break;
              case 'T':
                kmer_b.set(j * 2);
                kmer_b.set(j * 2 + 1);
                break;
              default:
                spdlog::error("unknown character: {}", kmer_s[j]);
                exit(1);
            }
          }

          int bucket;
          std::bitset<K * 2 - B> key;
          std::tie(bucket, key) = GetBucketAndKey(kmer_b);

          std::lock_guard lck{buckets_mutexes[bucket]};
          buckets_[bucket].insert(key);
        }
      }
    }
  }

  // Loads from a dumped data.
  KmerSet(std::istream& is) {
    int k, b;
    is >> k >> b;
    if (k != K || b != B) {
      spdlog::error("incompatible data");
      return;
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

  bool Contains(const std::bitset<K * 2>& kmer) const {
    int bucket;
    std::bitset<K * 2 - B> key;
    std::tie(bucket, key) = GetBucketAndKey(kmer);
    return buckets_[bucket].find(key) != buckets_[bucket].end();
  }

 private:
  std::array<absl::flat_hash_set<std::bitset<K * 2 - B>>, 1 << B> buckets_;

  static std::pair<int, std::bitset<K * 2 - B>> GetBucketAndKey(
      const std::bitset<K * 2>& kmer) {
    const int key_bits = K * 2 - B;
    const int bucket = kmer.to_ullong() >> key_bits;
    const std::bitset<K * 2 - B> key{kmer.to_ullong() % (1ull << key_bits)};
    return {bucket, key};
  }
};

std::string GetFormattedTime() {
  absl::Time now = absl::Now();
  return absl::FormatTime(now, absl::UTCTimeZone());
}

int main() {
  std::ios_base::sync_with_stdio(false);

  spdlog::info("constructing kmer_set");
  KmerSet kmer_set(std::cin, 100'000'000, 500'000'000);
  spdlog::info("constructed kmer_set");
  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());

  {
    spdlog::info("dumping kmer_set");
    std::ofstream kmer_dump_file;
    kmer_dump_file.open("kmer_dump.txt");
    kmer_set.Dump(kmer_dump_file);
    kmer_dump_file.close();
    spdlog::info("dumped kmer_set");
  }

  {
    spdlog::info("loading kmer_set_copy");
    std::ifstream kmer_dump_file{"kmer_dump.txt"};
    KmerSet kmer_set_copy(kmer_dump_file);
    spdlog::info("loaded kmer_set_copy");
    spdlog::info("kmer_set_copy.Size() = {}", kmer_set_copy.Size());
  }
}
