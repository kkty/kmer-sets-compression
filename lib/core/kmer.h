#ifndef CORE_KMER_H_
#define CORE_KMER_H_

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>

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
    std::uint64_t bits = 0;

    for (int i = 0; i < K; i++) {
      bits <<= 2;

      switch (s[i]) {
        case 'A':
          break;
        case 'C':
          bits += 1;
          break;
        case 'G':
          bits += 2;
          break;
        case 'T':
          bits += 3;
          break;
        default:
          assert(false);
      }
    }

    bits_ = bits;
  }

  explicit Kmer(std::uint64_t bits) : bits_(bits) {}

  // Returns the string representation of the kmer.
  // The length of the output string is K and only contains 'A', 'C', 'G' and
  // 'T'.
  std::string String() const {
    std::string s;
    s.resize(K);

    std::int64_t bits = bits_;

    for (int i = 0; i < K; i++) {
      switch (bits % 4) {
        case 0:
          s[K - 1 - i] = 'A';
          break;
        case 1:
          s[K - 1 - i] = 'C';
          break;
        case 2:
          s[K - 1 - i] = 'G';
          break;
        case 3:
          s[K - 1 - i] = 'T';
          break;
      }

      bits >>= 2;
    }

    return s;
  }

  char Last() const {
    switch (bits_ % 4) {
      case 0:
        return 'A';
      case 1:
        return 'C';
      case 2:
        return 'G';
      case 3:
        return 'T';
    }

    assert(false);
    return ' ';
  }

  // Returns the complement of the kmer.
  // The complement of a kmer is constructed by reversing the original kmer and
  // replacing 'A', 'C', 'G' and 'T' with 'T', 'G', 'C' and 'A' respectively.
  // Idea: kmer and kmer.Complement() can be paired in double stranded
  // structures.
  // Example: The complement of "AACCG" is "CGGTT".
  Kmer<K> Complement() const {
    Kmer<K> complement;

    std::int64_t bits = bits_;

    for (int i = 0; i < K; i++) {
      complement.bits_ <<= 2;

      switch (bits % 4) {
        case 0:
          complement.bits_ += 3;
          break;
        case 1:
          complement.bits_ += 2;
          break;
        case 2:
          complement.bits_ += 1;
          break;
        case 3:
          break;
      }

      bits >>= 2;
    }

    return complement;
  }

  // Returns the minimum of the kmer and the complement of the kmer.
  // The dictionary ordering is used to get the minimum.
  Kmer<K> Canonical() const { return std::min(*this, Complement()); }

  // Returns the concatenation of the (K-1)-suffix of the kmer and "c".
  Kmer<K> Next(char c) const {
    Kmer<K> next = *this;

    next.bits_ <<= 2;

    // Clears upper bits.
    next.bits_ &= std::numeric_limits<std::uint64_t>::max() >> (64 - K * 2);

    switch (c) {
      case 'A':
        break;
      case 'C':
        next.bits_ += 1;
        break;
      case 'G':
        next.bits_ += 2;
        break;
      case 'T':
        next.bits_ += 3;
        break;
      default:
        assert(false);
    }

    return next;
  }

  // Returns the concatenation of "c" and the (K-1)-prefix of the kmer.
  Kmer<K> Prev(char c) const {
    Kmer<K> prev = *this;

    prev.bits_ >>= 2;

    switch (c) {
      case 'A':
        break;
      case 'C':
        prev.bits_ += static_cast<std::uint64_t>(1) << ((K - 1) * 2);
        break;
      case 'G':
        prev.bits_ += static_cast<std::uint64_t>(2) << ((K - 1) * 2);
        break;
      case 'T':
        prev.bits_ += static_cast<std::uint64_t>(3) << ((K - 1) * 2);
        break;
      default:
        assert(false);
    }

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
  std::uint64_t Bits() const { return bits_; }

  std::size_t Hash() const { return bits_; }

  template <typename H>
  friend H AbslHashValue(H h, const Kmer& kmer) {
    return H::combine(std::move(h), kmer.bits_);
  }

 private:
  std::uint64_t bits_ = 0;
};

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
  return lhs.Bits() < rhs.Bits();
}

template <int K>
bool operator>(const Kmer<K>& lhs, const Kmer<K>& rhs) {
  return lhs.Bits() > rhs.Bits();
}

#endif