#include "core/kmer_set_immutable.h"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "core/kmer.h"
#include "gtest/gtest.h"

TEST(KmerSetImmutable, FromAndToKmers) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 100000;

  std::vector<Kmer<K>> kmers = GetRandomKmers<K>(n);

  KmerSetImmutable<K, N, KeyType> kmer_set_immutable(kmers, n_workers);

  std::vector<Kmer<K>> reconstructed = kmer_set_immutable.ToKmers(n_workers);

  std::sort(kmers.begin(), kmers.end());
  std::sort(reconstructed.begin(), reconstructed.end());

  ASSERT_EQ(kmers, reconstructed);
}
