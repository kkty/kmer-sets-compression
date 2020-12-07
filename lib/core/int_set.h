#ifndef CORE_INT_SET_H_
#define CORE_INT_SET_H_

#include <cstdlib>
#include <cstring>
#include <vector>

#include "streamvbyte.h"

// IntSet can be used to store an immutable set of integers.
// It is not possible to add or remove elements from the structure, but, fast
// operations are provided for adding two IntSets, intersecting two IntSets, and
// subtracting one IntSet from another.
template <typename T>
class IntSet {
 public:
  IntSet() = default;

  IntSet(const IntSet& other)
      : first_(other.first_),
        n_(other.n_),
        compressed_diff_(new uint8_t[other.compressed_diff_size_]),
        compressed_diff_size_(other.compressed_diff_size_) {
    std::memcpy(compressed_diff_, other.compressed_diff_,
                compressed_diff_size_);
  }

  IntSet(IntSet&& other) noexcept
      : first_(other.first_),
        n_(other.n_),
        compressed_diff_(other.compressed_diff_),
        compressed_diff_size_(other.compressed_diff_size_) {
    other.n_ = 0;
    other.compressed_diff_ = nullptr;
  }

  IntSet& operator=(const IntSet& other) {
    if (this == &other) return *this;

    delete[] compressed_diff_;

    first_ = other.first_;
    n_ = other.n_;
    compressed_diff_size_ = other.compressed_diff_size_;
    compressed_diff_ = new uint8_t[compressed_diff_size_];

    std::memcpy(compressed_diff_, other.compressed_diff_,
                compressed_diff_size_);

    return *this;
  }

  IntSet& operator=(IntSet&& other) noexcept {
    delete[] compressed_diff_;

    first_ = other.first_;
    n_ = other.n_;
    compressed_diff_ = other.compressed_diff_;
    compressed_diff_size_ = other.compressed_diff_size_;

    other.compressed_diff_ = nullptr;
    other.n_ = 0;

    return *this;
  }

  explicit IntSet(const std::vector<T>& v) : n_(v.size()) {
    if (n_ == 0) {
      return;
    }

    first_ = v.front();

    if (n_ == 1) {
      return;
    }

    std::vector<uint32_t> diff(n_ - 1);

    for (int64_t i = 0; i < n_ - 1; i++) {
      diff[i] = v[i + 1] - v[i];
    }

    compressed_diff_ = new uint8_t[streamvbyte_max_compressedbytes(n_ - 1)];

    compressed_diff_size_ =
        streamvbyte_encode(diff.data(), n_ - 1, compressed_diff_);
  }

  ~IntSet() {
    if (compressed_diff_ == nullptr) {
      return;
    }

    delete[] compressed_diff_;
  }

  // Reconstructs the original vector.
  std::vector<T> Decode() const {
    std::vector<T> v;
    ForEach([&](T i) { v.push_back(i); });
    return v;
  }

  // Returns the number of elements in the set.
  int64_t Size() const { return n_; }

  // Executes "func" for each element in the set.
  template <typename FuncType>
  void ForEach(FuncType func) const {
    if (n_ == 0) return;

    T current = first_;
    std::vector<uint32_t> diff = Diff();
    auto it = diff.begin();

    while (true) {
      func(current);
      if (it == diff.end()) return;
      current += *it;
      it++;
    }
  }

  // Executes "func" for each element that is present in both of the two
  // IntSets.
  template <typename FuncType>
  void Intersection(const IntSet& other, FuncType func) const {
    const IntSet& lhs = *this;
    const IntSet& rhs = other;

    if (lhs.n_ == 0 || rhs.n_ == 0) return;

    std::vector<uint32_t> diff_lhs = lhs.Diff();
    std::vector<uint32_t> diff_rhs = rhs.Diff();

    auto it_lhs = diff_lhs.begin();
    auto it_rhs = diff_rhs.begin();

    T current_lhs = lhs.first_;
    T current_rhs = rhs.first_;

    while (true) {
      if (current_lhs < current_rhs) {
        if (it_lhs == diff_lhs.end()) return;

        current_lhs += *it_lhs;
        it_lhs++;
      } else if (current_lhs > current_rhs) {
        if (it_rhs == diff_rhs.end()) return;

        current_rhs += *it_rhs;
        it_rhs++;
      } else {
        func(current_lhs);

        if (it_lhs == diff_lhs.end() || it_rhs == diff_rhs.end()) return;

        current_lhs += *it_lhs;
        it_lhs++;

        current_rhs += *it_rhs;
        it_rhs++;
      }
    }
  }

  // Returns the size of the overlap between the two IntSets.
  int64_t IntersectionSize(const IntSet& other) const {
    int64_t size = 0;
    Intersection(other, [&](T) { size++; });
    return size;
  }

  // Returns the overlap between the two IntSets.
  IntSet Intersection(const IntSet& other) const {
    std::vector<T> v;
    Intersection(other, [&](T i) { v.push_back(i); });
    return IntSet(v);
  }

  // Executes "func" for each element that is present in either of the two
  // IntSets.
  template <typename FuncType>
  void Add(const IntSet& other, FuncType func) const {
    const IntSet& lhs = *this;
    const IntSet& rhs = other;

    if (lhs.n_ == 0 && rhs.n_ == 0) return;

    if (lhs.n_ == 0) {
      rhs.ForEach(func);
      return;
    }

    if (rhs.n_ == 0) {
      lhs.ForEach(func);
      return;
    }

    std::vector<uint32_t> diff_lhs = lhs.Diff();
    std::vector<uint32_t> diff_rhs = rhs.Diff();

    auto it_lhs = diff_lhs.begin();
    auto it_rhs = diff_rhs.begin();

    T current_lhs = lhs.first_;
    T current_rhs = rhs.first_;

    while (true) {
      if (current_lhs < current_rhs) {
        func(current_lhs);

        if (it_lhs == diff_lhs.end()) {
          while (true) {
            func(current_rhs);

            if (it_rhs == diff_rhs.end()) return;

            current_rhs += *it_rhs;
            it_rhs++;
          }
        }

        current_lhs += *it_lhs;
        it_lhs++;
      } else if (current_lhs > current_rhs) {
        func(current_rhs);

        if (it_rhs == diff_rhs.end()) {
          while (true) {
            func(current_lhs);

            if (it_lhs == diff_lhs.end()) return;

            current_lhs += *it_lhs;
            it_lhs++;
          }
        }

        current_rhs += *it_rhs;
        it_rhs++;
      } else {
        func(current_lhs);

        if (it_lhs == diff_lhs.end() && it_rhs == diff_rhs.end()) return;

        if (it_lhs == diff_lhs.end()) {
          while (true) {
            if (it_rhs == diff_rhs.end()) return;

            current_rhs += *it_rhs;
            it_rhs++;

            func(current_rhs);
          }
        }

        if (it_rhs == diff_rhs.end()) {
          while (true) {
            if (it_lhs == diff_lhs.end()) return;

            current_lhs += *it_lhs;
            it_lhs++;

            func(current_lhs);
          }
        }

        current_lhs += *it_lhs;
        it_lhs++;

        current_rhs += *it_rhs;
        it_rhs++;
      }
    }
  }

  // Adds two IntSets.
  IntSet Add(const IntSet& other) const {
    std::vector<T> v;
    Add(other, [&](T i) { v.push_back(i); });
    return IntSet(v);
  }

  // Executes "func" for each element that is present in one IntSet, but not in
  // the other.
  template <typename FuncType>
  void Sub(const IntSet& other, FuncType func) const {
    const IntSet& lhs = *this;
    const IntSet& rhs = other;

    if (lhs.n_ == 0) return;

    if (rhs.n_ == 0) {
      lhs.ForEach(func);
      return;
    }

    std::vector<uint32_t> diff_lhs = lhs.Diff();
    std::vector<uint32_t> diff_rhs = rhs.Diff();

    auto it_lhs = diff_lhs.begin();
    auto it_rhs = diff_rhs.begin();

    T current_lhs = lhs.first_;
    T current_rhs = rhs.first_;

    while (true) {
      if (current_lhs < current_rhs) {
        func(current_lhs);

        if (it_lhs == diff_lhs.end()) return;

        current_lhs += *it_lhs;
        it_lhs++;
      } else if (current_lhs > current_rhs) {
        if (it_rhs == diff_rhs.end()) {
          while (true) {
            func(current_lhs);

            if (it_lhs == diff_lhs.end()) return;

            current_lhs += *it_lhs;
            it_lhs++;
          }
        }

        current_rhs += *it_rhs;
        it_rhs++;
      } else {
        if (it_lhs == diff_lhs.end()) return;

        if (it_rhs == diff_rhs.end()) {
          while (true) {
            if (it_lhs == diff_lhs.end()) return;

            current_lhs += *it_lhs;
            it_lhs++;

            func(current_lhs);
          }
        }

        current_lhs += *it_lhs;
        it_lhs++;

        current_rhs += *it_rhs;
        it_rhs++;
      }
    }
  }

  // Subtracts between two IntSets.
  IntSet Sub(const IntSet& other) const {
    std::vector<T> v;
    Sub(other, [&](T i) { v.push_back(i); });
    return IntSet(v);
  }

 private:
  T first_;
  int64_t n_ = 0;

  std::vector<uint32_t> Diff() const {
    assert(n_ > 0);
    std::vector<uint32_t> diff(n_ - 1);
    streamvbyte_decode(compressed_diff_, diff.data(), n_ - 1);
    return diff;
  }

  uint8_t* compressed_diff_ = nullptr;
  int64_t compressed_diff_size_ = 0;
};

#endif
