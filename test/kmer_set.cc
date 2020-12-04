#include "core/kmer_set.h"

#include <bitset>
#include <tuple>
#include <vector>

#include "core/kmer.h"
#include "core/kmer_counter.h"
#include "gtest/gtest.h"

TEST(KmerSet, BucketAndKey) {
  const int K = 5;
  using KeyType = uint8_t;

  Kmer<K> kmer("AGCTG");

  int bucket;
  KeyType key;
  std::tie(bucket, key) = GetBucketAndKeyFromKmer<K, KeyType>(kmer);

  ASSERT_EQ(kmer.String(),
            (GetKmerFromBucketAndKey<K, KeyType>(bucket, key)).String());
}

TEST(KmerSet, AddRemove) {
  const int K = 5;
  using KeyType = uint8_t;

  Kmer<K> kmer("AAAAA");

  KmerSet<K, KeyType> kmer_set;

  ASSERT_EQ(kmer_set.Size(), 0);
  ASSERT_FALSE(kmer_set.Contains(kmer));

  kmer_set.Add(kmer);

  ASSERT_EQ(kmer_set.Size(), 1);
  ASSERT_TRUE(kmer_set.Contains(kmer));

  kmer_set.Remove(kmer);

  ASSERT_EQ(kmer_set.Size(), 0);
  ASSERT_FALSE(kmer_set.Contains(kmer));
}

TEST(KmerSet, Find) {
  const int K = 5;
  using KeyType = uint8_t;

  KmerSet<K, KeyType> kmer_set;
  kmer_set.Add(Kmer<K>("AAAAA"));
  kmer_set.Add(Kmer<K>("CCCCC"));

  {
    const std::vector<Kmer<K>> kmers = kmer_set.Find(
        [](const Kmer<K>& kmer) { return kmer.String()[0] == 'A'; }, 1);
    ASSERT_EQ(kmers.size(), 1);
    ASSERT_EQ(kmers[0].String(), "AAAAA");
  }

  {
    const std::vector<Kmer<K>> kmers = kmer_set.Find(
        [](const Kmer<K>& kmer) { return kmer.String()[1] == 'C'; }, 1);
    ASSERT_EQ(kmers.size(), 1);
    ASSERT_EQ(kmers[0].String(), "CCCCC");
  }
}

TEST(KmerSet, Operators) {
  const int K = 5;
  using KeyType = uint8_t;

  KmerSet<K, KeyType> kmer_set1;
  KmerSet<K, KeyType> kmer_set2;

  kmer_set1.Add(Kmer<K>("AAAAA"));
  kmer_set1.Add(Kmer<K>("TTTTT"));
  kmer_set1.Add(Kmer<K>("CCCCC"));

  kmer_set2.Add(Kmer<K>("AAAAA"));
  kmer_set2.Add(Kmer<K>("TTTTT"));
  kmer_set2.Add(Kmer<K>("GGGGG"));

  ASSERT_EQ(Add(kmer_set1, kmer_set2, 1).Size(), 4);
  ASSERT_EQ(Sub(kmer_set1, kmer_set2, 1).Size(), 1);
  ASSERT_EQ(Sub(kmer_set2, kmer_set1, 1).Size(), 1);
  ASSERT_EQ(Sub(kmer_set2, kmer_set1, 1).Size(), 1);
  ASSERT_EQ(Intersection(kmer_set2, kmer_set1, 1).Size(), 2);
  ASSERT_EQ(Intersection(kmer_set2, kmer_set1, 1).Size(), 2);
}

TEST(KmerSet, EqualsAndDiff) {
  const int K = 5;
  using KeyType = uint8_t;

  KmerSet<K, KeyType> kmer_set1;
  KmerSet<K, KeyType> kmer_set2;
  KmerSet<K, KeyType> kmer_set3;

  kmer_set1.Add(Kmer<K>("AAAAA"));
  kmer_set1.Add(Kmer<K>("TTTTT"));

  kmer_set2.Add(Kmer<K>("AAAAA"));
  kmer_set2.Add(Kmer<K>("TTTTT"));

  kmer_set3.Add(Kmer<K>("AAAAA"));
  kmer_set3.Add(Kmer<K>("CCCCC"));
  kmer_set3.Add(Kmer<K>("GGGGG"));

  ASSERT_TRUE(kmer_set1.Equals(kmer_set1, 1));
  ASSERT_TRUE(kmer_set2.Equals(kmer_set2, 1));
  ASSERT_TRUE(kmer_set3.Equals(kmer_set3, 1));

  ASSERT_TRUE(kmer_set1.Equals(kmer_set2, 1));
  ASSERT_TRUE(kmer_set2.Equals(kmer_set1, 1));

  ASSERT_FALSE(kmer_set1.Equals(kmer_set3, 1));
  ASSERT_FALSE(kmer_set3.Equals(kmer_set1, 1));

  ASSERT_EQ(kmer_set1.Diff(kmer_set3, 1), 3);
  ASSERT_EQ(kmer_set3.Diff(kmer_set1, 1), 3);

  ASSERT_EQ(kmer_set1.Common(kmer_set3, 1), 1);
  ASSERT_EQ(kmer_set3.Common(kmer_set1, 1), 1);

  ASSERT_DOUBLE_EQ(kmer_set1.Similarity(kmer_set3, 1), 0.25);
  ASSERT_DOUBLE_EQ(kmer_set3.Similarity(kmer_set1, 1), 0.25);
}

TEST(KmerSet, Extract) {
  const int K = 5;
  using KeyType = uint8_t;

  KmerSet<K, KeyType> kmer_set;
  kmer_set.Add(Kmer<K>("ACGTA"));
  kmer_set.Add(Kmer<K>("CGTAC"));

  const int K2 = 4;
  using KeyType2 = uint8_t;

  KmerSet<K2, KeyType2> extracted = kmer_set.template Extract<K2, KeyType2>(1);
  ASSERT_EQ(extracted.Size(), 3);
  ASSERT_TRUE(extracted.Contains(Kmer<K2>("ACGT")));
  ASSERT_TRUE(extracted.Contains(Kmer<K2>("CGTA")));
  ASSERT_TRUE(extracted.Contains(Kmer<K2>("GTAC")));
}
