#ifndef CORE_KMER_H_
#define CORE_KMER_H_

#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"

// Kmer is used to represent a single kmer.
// Data is represented with 2 * K bits internally, regarding 'A', 'C', 'G', and
// 'T' as 00, 01, 10 and 11 respectively.
// The class supports conversions from / to strings, which may contain 'A', 'C',
// 'G' and 'T'.
template <int K>
class Kmer {
 public:
  Kmer() = default;

  explicit Kmer(const std::string& s) {
    for (int i = 0; i < K; i++) {
      assert(s[i] == 'A' || s[i] == 'C' || s[i] == 'G' || s[i] == 'T');
      Set(i, s[i]);
    }
  }

  explicit Kmer(const std::bitset<K * 2>& bits) : bits_(bits) {}

  // Returns the string representation of the kmer.
  // The length of the output string is K and only contains 'A', 'C', 'G' and
  // 'T'.
  std::string String() const {
    std::string s;
    s.reserve(K);
    for (int i = 0; i < K; i++) {
      s += Get(i);
    }
    return s;
  }

  // Returns the ith character of the kmer. The return value is either 'A', 'C',
  // 'G' or 'T'.
  char Get(int idx) const {
    // Gets the ith bit from the right.
    const auto get = [&](int i) { return bits_[K * 2 - 1 - i]; };

    const bool b0 = get(idx * 2);
    const bool b1 = get(idx * 2 + 1);

    if (!b0) {
      if (!b1) {
        // 00
        return 'A';
      } else {
        // 01
        return 'C';
      }
    } else {
      if (!b1) {
        // 10
        return 'G';
      } else {
        // 11
        return 'T';
      }
    }
  }

  // Sets the ith character of the kmer. "idx" should be in [0, K) and "c"
  // should be either 'A', 'C', 'G' or 'T'.
  void Set(int idx, char c) {
    // Sets the ith bit from the left.
    const auto set = [&](int i) { bits_.set(K * 2 - 1 - i); };

    switch (c) {
      case 'A':
        break;
      case 'C':
        set(idx * 2 + 1);
        break;
      case 'G':
        set(idx * 2);
        break;
      case 'T':
        set(idx * 2);
        set(idx * 2 + 1);
        break;
      default:
        assert(false);
    }
  }

  // Returns the complement of the kmer.
  // The complement of a kmer is constructed by reversing the original kmer and
  // replacing 'A', 'C', 'G' and 'T' with 'T', 'G', 'C' and 'A' respectively.
  // Idea: kmer and kmer.Complement() can be paired in double stranded
  // structures.
  // Example: The complement of "AACCG" is "CGGTT".
  Kmer<K> Complement() const {
    Kmer<K> complement;
    for (int i = 0; i < K; i++) {
      complement.Set(K - 1 - i, Get(i));
    }
    complement.bits_.flip();
    return complement;
  }

  // Returns the minimum of the kmer and the complement of the kmer.
  // The dictionary ordering is used to get the minimum.
  Kmer<K> Canonical() const { return std::min(*this, Complement()); }

  // Returns the concatenation of the (K-1)-suffix of the kmer and "c".
  Kmer<K> Next(char c) const {
    Kmer<K> next = *this;
    next.bits_ <<= 2;
    next.Set(K - 1, c);
    return next;
  }

  // Returns the concatenation of "c" and the (K-1)-prefix of the kmer.
  Kmer<K> Prev(char c) const {
    Kmer<K> prev = *this;
    prev.bits_ >>= 2;
    prev.Set(0, c);
    return prev;
  }

  // Returns {Next('A'), Next('C'), Next('G'), Next('T')}.
  std::array<Kmer<K>, 4> Nexts() const {
    std::array<Kmer<K>, 4> nexts;
    nexts[0] = Next('A');
    nexts[1] = Next('C');
    nexts[2] = Next('G');
    nexts[3] = Next('T');
    return nexts;
  }

  // Returns {Prev('A'), Prev('C'), Prev('G'), Prev('T')}.
  std::array<Kmer<K>, 4> Prevs() const {
    std::array<Kmer<K>, 4> prevs;
    prevs[0] = Prev('A');
    prevs[1] = Prev('C');
    prevs[2] = Prev('G');
    prevs[3] = Prev('T');
    return prevs;
  }

  // Returns the bit-wise representation of the kmer.
  std::bitset<K * 2> Bits() const { return bits_; }

  std::size_t Hash() const { return std::hash<std::bitset<K * 2>>()(bits_); }

 private:
  std::bitset<K * 2> bits_;
};

// Returns a randomly-generated kmer.
template <int K>
Kmer<K> GetRandomKmer() {
  absl::InsecureBitGen bitgen;
  Kmer<K> kmer;

  for (int i = 0; i < K; i++) {
    switch (absl::Uniform(bitgen, 0, 4)) {
      case 0:
        kmer.Set(i, 'A');
        break;
      case 1:
        kmer.Set(i, 'G');
        break;
      case 2:
        kmer.Set(i, 'C');
        break;
      case 3:
        kmer.Set(i, 'T');
        break;
    }
  }

  return kmer;
}

// Returns n randomly-generated kmers.
template <int K>
std::vector<Kmer<K>> GetRandomKmers(int n) {
  absl::flat_hash_set<Kmer<K>> kmers;

  while (static_cast<int>(kmers.size()) < n) {
    kmers.insert(GetRandomKmer<K>());
  }

  return std::vector<Kmer<K>>(kmers.begin(), kmers.end());
}

template <int K>
bool operator==(const Kmer<K>& lhs, const Kmer<K>& rhs) {
  return lhs.Bits() == rhs.Bits();
}

template <int K>
bool operator!=(const Kmer<K>& lhs, const Kmer<K>& rhs) {
  return !(lhs == rhs);
}

template <int K>
bool operator<(const Kmer<K>& lhs, const Kmer<K>& rhs) {
  return lhs.Bits().to_ullong() < rhs.Bits().to_ullong();
}

template <int K>
bool operator>(const Kmer<K>& lhs, const Kmer<K>& rhs) {
  return lhs.Bits().to_ullong() > rhs.Bits().to_ullong();
}

namespace std {
template <int K>
struct hash<Kmer<K>> {
  std::size_t operator()(const Kmer<K>& kmer) const { return kmer.Hash(); }
};
}  // namespace std

#endif