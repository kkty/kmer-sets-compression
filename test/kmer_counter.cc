#include "kmer_counter.h"

#include <tuple>

#include "gtest/gtest.h"
#include "kmer.h"

TEST(KmerCounter, AddAndGet) {
  const int K = 5;
  using KeyType = uint8_t;
  using ValueType = uint8_t;

  KmerCounter<K, KeyType, ValueType> kmer_counter;

  const Kmer<K> kmer1("AAAAA");
  const Kmer<K> kmer2("CCCCC");
  const Kmer<K> kmer3("TTTTT");

  kmer_counter.Add(kmer1, 1);
  kmer_counter.Add(kmer2, 2);
  kmer_counter.Add(kmer3, 3);
  kmer_counter.Add(kmer1, 1);

  ASSERT_EQ(kmer_counter.Get(kmer1), 2);
  ASSERT_EQ(kmer_counter.Get(kmer2), 2);
  ASSERT_EQ(kmer_counter.Get(kmer3), 3);
}

TEST(KmerCounter, Multiply) {
  const int K = 5;
  const int n_workers = 1;
  using KeyType = uint8_t;
  using ValueType = uint8_t;

  KmerCounter<K, KeyType, ValueType> kmer_counter;

  const Kmer<K> kmer1("AAAAA");
  const Kmer<K> kmer2("CCCCC");

  kmer_counter.Add(kmer1, 1);
  kmer_counter.Add(kmer2, 2);

  kmer_counter.Multiply(3, n_workers);

  ASSERT_EQ(kmer_counter.Get(kmer1), 3);
  ASSERT_EQ(kmer_counter.Get(kmer2), 6);
}

TEST(KmerCounter, ToSetAndFromSet) {
  const int K = 5;
  const int n_workers = 1;
  using KeyType = uint8_t;
  using ValueType = uint8_t;

  KmerCounter<K, KeyType, ValueType> kmer_counter;

  kmer_counter.Add(Kmer<K>("AAAAA"), 3);
  kmer_counter.Add(Kmer<K>("CCCCC"), 1);
  kmer_counter.Add(Kmer<K>("GGGGG"), 2);
  kmer_counter.Add(Kmer<K>("TTTTT"), 4);

  KmerSet<K, KeyType> kmer_set;
  int cutoff;
  std::tie(kmer_set, cutoff) = kmer_counter.ToSet(3, n_workers);

  ASSERT_EQ(cutoff, 2);
  ASSERT_EQ(kmer_set.Size(), 2);
  ASSERT_TRUE(kmer_set.Contains(Kmer<K>("AAAAA")));
  ASSERT_TRUE(kmer_set.Contains(Kmer<K>("TTTTT")));

  KmerCounter<K, KeyType, ValueType> reconstructed_kmer_counter =
      KmerCounter<K, KeyType, ValueType>::FromSet(kmer_set, n_workers);

  ASSERT_EQ(reconstructed_kmer_counter.Size(), 2);
  ASSERT_EQ(reconstructed_kmer_counter.Get(Kmer<K>("AAAAA")), n_workers);
  ASSERT_EQ(reconstructed_kmer_counter.Get(Kmer<K>("TTTTT")), n_workers);
}