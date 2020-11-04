#ifndef KMER_H_
#define KMER_H_

#include <array>
#include <bitset>
#include <string>

template <int K>
class Kmer {
 public:
  Kmer() = default;

  Kmer(const std::string& s) {
    for (int i = 0; i < K; i++) Set(i, s[i]);
  }

  Kmer(const std::bitset<K * 2>& bits) : bits_(bits) {}

  std::string String() const {
    std::string s;
    s.reserve(K);
    for (int i = 0; i < K; i++) {
      s += Get(i);
    }
    return s;
  }

  char Get(int idx) const {
    // Gets the ith bit from the right.
    const auto get = [&](int i) { return bits_[K * 2 - 1 - i]; };

    const bool b0 = get(idx * 2);
    const bool b1 = get(idx * 2 + 1);

    if (!b0 && !b1) return 'A';
    if (b0 && !b1) return 'C';
    if (!b0 && b1) return 'G';
    return 'T';
  }

  void Set(int idx, char c) {
    // Sets the ith bit from the left.
    const auto set = [&](int i) { bits_.set(K * 2 - 1 - i); };

    switch (c) {
      case 'A':
        break;
      case 'C':
        set(idx * 2);
        break;
      case 'G':
        set(idx * 2 + 1);
        break;
      case 'T':
        set(idx * 2);
        set(idx * 2 + 1);
        break;
    }
  }

  Kmer<K> Complement() const {
    Kmer<K> complement = *this;
    complement.bits_.flip();
    return complement;
  }

  Kmer<K> Next(char c) const {
    Kmer<K> next = *this;
    next.bits_ <<= 2;
    next.Set(K - 1, c);
    return next;
  }

  Kmer<K> Prev(char c) const {
    Kmer<K> prev = *this;
    prev.bits_ >>= 2;
    prev.Set(0, c);
    return prev;
  }

  std::array<Kmer<K>, 4> Nexts() const {
    std::array<Kmer<K>, 4> nexts;
    nexts[0] = Next('A');
    nexts[1] = Next('C');
    nexts[2] = Next('G');
    nexts[3] = Next('T');
    return nexts;
  }

  std::array<Kmer<K>, 4> Prevs() const {
    std::array<Kmer<K>, 4> prevs;
    prevs[0] = Prev('A');
    prevs[1] = Prev('C');
    prevs[2] = Prev('G');
    prevs[3] = Prev('T');
    return prevs;
  }

  std::bitset<K * 2> Bits() const { return bits_; }

  size_t Hash() const { return std::hash<std::bitset<K * 2>>()(bits_); }

 private:
  std::bitset<K * 2> bits_;
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
  return lhs.Bits().to_ullong() < rhs.Bits().to_ullong();
}

namespace std {
template <int K>
struct hash<Kmer<K>> {
  size_t operator()(const Kmer<K>& kmer) const { return kmer.Hash(); }
};
}  // namespace std

#endif