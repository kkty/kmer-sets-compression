#ifndef CORE_RANDOM_H_
#define CORE_RANDOM_H_

#include <algorithm>
#include <limits>
#include <vector>

#include "absl/container/flat_hash_set.h"
#include "absl/random/random.h"

// Returns n randomly generated integers from min to max.
template <typename T = int>
std::vector<T> GetRandomInts(int n, bool is_unique = false,
                             bool is_sorted = false,
                             T min = std::numeric_limits<T>::min(),
                             T max = std::numeric_limits<T>::max()) {
  absl::InsecureBitGen bitgen;

  std::vector<T> v;
  v.reserve(n);

  if (is_unique) {
    absl::flat_hash_set<T> s;

    while (static_cast<int>(s.size()) < n) {
      s.insert(absl::Uniform(absl::IntervalClosed, bitgen, min, max));
    }

    v.insert(v.end(), s.begin(), s.end());
  } else {
    for (int i = 0; i < n; i++) {
      v.push_back(absl::Uniform(absl::IntervalClosed, bitgen, min, max));
    }
  }

  if (is_sorted) {
    std::sort(v.begin(), v.end());
  }

  return v;
}

#endif