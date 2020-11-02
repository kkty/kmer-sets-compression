#include "kmer_set.h"

#include <bitset>
#include <tuple>

#include "gtest/gtest.h"
#include "kmer.h"
#include "kmer_counter.h"

TEST(kmer_set, BucketAndKey) {
  const int K = 5;
  const int N = 2;
  Kmer<K> kmer("AGCTG");
  const auto [bucket, key] = GetBucketAndKeyFromKmer<K, N>(kmer);
  ASSERT_EQ(kmer, (GetKmerFromBucketAndKey<K, N>(bucket, key)));
}
