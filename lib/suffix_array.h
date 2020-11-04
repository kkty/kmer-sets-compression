#ifndef SUFFIX_ARRAY_H_
#define SUFFIX_ARRAY_H_

#include <algorithm>
#include <string>
#include <vector>

class SuffixArray {
 public:
  SuffixArray() = default;

  SuffixArray(const std::string& s) : s_(s), v_(s.length()) {
    for (int64_t i = 0; i < (int64_t)s_.length(); i++) v_[i] = i;

    std::sort(v_.begin(), v_.end(), [&](int64_t x, int64_t y) {
      // Equivalent to:
      // return s_.substr(x, s_.length() - x) < s_.substr(y, s_.length() - y);

      for (int64_t i = 0; i + std::max(x, y) < (int64_t)s_.length(); i++) {
        if (s_[x + i] < s_[y + i]) return true;
        if (s_[x + i] > s_[y + i]) return false;
      }

      return x > y;
    });
  }

  // Returns a list of indexes in the original string.
  std::vector<int64_t> Find(const std::string& query) const {
    // begin == -1 or s_[v_[begin]:] < query.
    int begin = -1;

    // end == s_.length() or s_[v_[end]:] >= query.
    int end = s_.length();

    while (end - begin > 1) {
      const int mid = (begin + end) / 2;

      // Comparison between s_[v_[mid]:] and s.
      // "less" is true if s_[v_[mid]:] < query.
      const bool less = [&] {
        for (int64_t i = 0;
             i + v_[mid] < (int64_t)s_.length() && i < (int64_t)query.length();
             i++) {
          const char c_s = s_[i + v_[mid]];
          const char c_query = query[i];
          if (c_s < c_query) return true;
          if (c_s > c_query) return false;
        }
        return s_.length() - v_[mid] < query.length();
      }();

      if (less)
        begin = mid;
      else
        end = mid;
    }

    // Returns true if query is a prefix of s_[v_[i]].
    const auto is_prefix = [&](int64_t i) {
      if ((int64_t)s_.length() - v_[i] < (int64_t)query.length()) return false;
      for (int j = 0; j < (int)query.length(); j++) {
        if (s_[v_[i] + j] != query[j]) return false;
      }
      return true;
    };

    std::vector<int64_t> indices;

    for (int64_t i = end; i < (int64_t)s_.length() && is_prefix(i); i++) {
      indices.push_back(v_[i]);
    }

    std::sort(indices.begin(), indices.end());

    return indices;
  }

  // Returns the original string.
  const std::string& String() const { return s_; }

 private:
  std::string s_;
  std::vector<int64_t> v_;
};

#endif