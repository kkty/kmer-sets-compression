#include "core/kmer_counter.h"

#include <cstdint>
#include <limits>
#include <tuple>
#include <utility>

#include "core/kmer.h"
#include "core/kmer_set.h"
#include "gtest/gtest.h"

TEST(KmerCounter, AddWithMax) {
  std::uint8_t x = 1;
  x = AddWithMax(std::numeric_limits<std::uint8_t>::max(), x);
  ASSERT_EQ(x, std::numeric_limits<std::uint8_t>::max());
}

TEST(KmerCounter, AddAndGet) {
  const int K = 5;
  const int N = 3;
  using KeyType = std::uint8_t;
  using ValueType = std::uint8_t;

  KmerCounter<K, N, KeyType, ValueType> kmer_counter;

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

TEST(KmerCounter, ToKmerSet) {
  const int K = 5;
  const int N = 3;
  const int n_workers = 1;
  using KeyType = std::uint8_t;
  using ValueType = std::uint8_t;

  KmerCounter<K, N, KeyType, ValueType> kmer_counter;

  kmer_counter.Add(Kmer<K>("AAAAA"), 3);
  kmer_counter.Add(Kmer<K>("CCCCC"), 1);
  kmer_counter.Add(Kmer<K>("GGGGG"), 2);
  kmer_counter.Add(Kmer<K>("TTTTT"), 4);

  KmerSet<K, N, KeyType> kmer_set;
  int cutoff;
  std::tie(kmer_set, cutoff) = kmer_counter.ToKmerSet(3, n_workers);

  ASSERT_EQ(cutoff, 2);
  ASSERT_EQ(kmer_set.Size(), 2);
  ASSERT_TRUE(kmer_set.Contains(Kmer<K>("AAAAA")));
  ASSERT_TRUE(kmer_set.Contains(Kmer<K>("TTTTT")));
}

TEST(KmerCounter, FromReads) {
  const int K = 5;
  const int N = 3;
  using KeyType = std::uint8_t;
  using ValueType = std::uint8_t;
  const int n_workers = 1;

  const std::vector<std::string> reads{
      "AACCGTT",
      "AACCGTA",
  };

  KmerCounter<K, N, KeyType, ValueType> kmer_counter;

  {
    absl::StatusOr<KmerCounter<K, N, KeyType, ValueType>> statusor =
        KmerCounter<K, N, KeyType, ValueType>::FromReads(reads, false,
                                                         n_workers);

    ASSERT_TRUE(statusor.ok());
    kmer_counter = std::move(statusor).value();
  }

  ASSERT_EQ(kmer_counter.Get(Kmer<K>("AACCG")), 2);
  ASSERT_EQ(kmer_counter.Get(Kmer<K>("ACCGT")), 2);
  ASSERT_EQ(kmer_counter.Get(Kmer<K>("CCGTT")), 1);
  ASSERT_EQ(kmer_counter.Get(Kmer<K>("CCGTA")), 1);
}
