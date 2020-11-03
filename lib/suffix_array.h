#ifndef SUFFIX_ARRAY_H_
#define SUFFIX_ARRAY_H_

#include <algorithm>
#include <string>
#include <vector>

class SuffixArray {
 public:
  SuffixArray() = default;

  SuffixArray(const std::string& s) : s_(s), v_(s.length()) {
    for (int i = 0; i < s_.length(); i++) v_[i] = i;

    std::sort(v_.begin(), v_.end(), [&](int x, int y) {
      // Equivalent to:
      // return s_.substr(x, s_.length() - x) < s_.substr(y, s_.length() - y);

      for (int i = 0; i + std::max(x, y) < s_.length(); i++) {
        if (s_[x + i] < s_[y + i]) return true;
        if (s_[x + i] > s_[y + i]) return false;
      }

      return x > y;
    });
  }

  // Returns a list of indexes in the original string.
  std::vector<int> Find(const std::string& query) const {
    // begin == -1 or s_[v_[begin]:] < query.
    int begin = -1;

    // end == s_.length() or s_[v_[end]:] >= query.
    int end = s_.length();

    while (end - begin > 1) {
      const int mid = (begin + end) / 2;

      // Comparison between s_[v_[mid]:] and s.
      // "less" is true if s_[v_[mid]:] < query.
      const bool less = [&] {
        for (int i = 0; i + v_[mid] < s_.length() && i < query.length(); i++) {
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
    const auto is_prefix = [&](int i) {
      if (s_.length() - v_[i] < query.length()) return false;
      for (int j = 0; j < (int)query.length(); j++) {
        if (s_[v_[i] + j] != query[j]) return false;
      }
      return true;
    };

    std::vector<int> indices;

    for (int i = end; i < s_.length() && is_prefix(i); i++) {
      indices.push_back(v_[i]);
    }

    std::sort(indices.begin(), indices.end());

    return indices;
  }

 private:
  std::string s_;
  std::vector<int> v_;
};

#endif