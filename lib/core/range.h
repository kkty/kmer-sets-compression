#ifndef CORE_RANGE_H_
#define CORE_RANGE_H_

#include <cstdint>
#include <vector>

class Range {
 public:
  class Iterator {
   public:
    explicit Iterator(std::int64_t current) : current_(current) {}

    std::int64_t operator*() const { return current_; }

    Iterator& operator++() {
      current_ += 1;
      return *this;
    }

    bool operator!=(const Iterator& other) const {
      return current_ != other.current_;
    }

   private:
    std::int64_t current_;
  };

  Range(std::int64_t begin, std::int64_t end) : begin_(begin), end_(end) {}

  Iterator begin() const { return Iterator(begin_); }

  Iterator end() const { return Iterator(end_); }

  std::int64_t Begin() const { return begin_; }

  std::int64_t End() const { return end_; }

  std::int64_t Size() const { return end_ - begin_; }

  // Splits the range to n distinct ranges. If n is smaller than Size(), some
  // resulting ranges have a size of zero.
  std::vector<Range> Split(std::int64_t n) const {
    const std::int64_t small_chunk_size = Size() / n;
    const std::int64_t large_chunk_size = small_chunk_size + 1;

    const std::int64_t large_chunk_n = Size() - small_chunk_size * n;
    const std::int64_t small_chunk_n = n - large_chunk_n;

    std::vector<Range> ranges;
    ranges.reserve(n);

    for (std::int64_t i = 0; i < small_chunk_n; i++) {
      if (i == 0) {
        ranges.emplace_back(begin_, begin_ + small_chunk_size);
      } else {
        ranges.emplace_back(ranges.back().end_,
                            ranges.back().end_ + small_chunk_size);
      }
    }

    for (std::int64_t i = 0; i < large_chunk_n; i++) {
      ranges.emplace_back(ranges.back().end_,
                          ranges.back().end_ + large_chunk_size);
    }

    return ranges;
  }

  // Executes a function for each integer in the range.
  // Range(...).ForEach([&](std::int64_t i) { ... });
  template <typename F>
  void ForEach(F f) const {
    for (std::int64_t i = begin_; i < end_; i++) {
      f(i);
    }
  }

 private:
  const std::int64_t begin_;
  const std::int64_t end_;
};

#endif