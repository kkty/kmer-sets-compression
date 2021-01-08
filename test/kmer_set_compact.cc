#include "core/kmer_set_compact.h"

#include <cstdint>
#include <utility>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "core/kmer_set.h"
#include "gtest/gtest.h"
#include "io.h"
#include "random.h"

TEST(KmerSetCompact, DumpAndLoad) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 100000;

  const KmerSet<K, N, KeyType> kmer_set =
      GetRandomKmerSet<K, N, KeyType>(n, true);

  const KmerSetCompact<K, N, KeyType> kmer_set_compact =
      KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, true, true,
                                                 n_workers);

  const TemporaryFile temporary_file;

  {
    const absl::Status status =
        kmer_set_compact.Dump(temporary_file.Name(), "", n_workers);

    ASSERT_TRUE(status.ok());
  }

  KmerSet<K, N, KeyType> reconstructed;

  {
    absl::StatusOr<KmerSetCompact<K, N, KeyType>> statusor =
        KmerSetCompact<K, N, KeyType>::Load(temporary_file.Name(), "",
                                            n_workers);

    ASSERT_TRUE(statusor.ok());

    reconstructed = std::move(statusor).value().ToKmerSet(true, n_workers);
  }

  ASSERT_TRUE(kmer_set.Equals(reconstructed, n_workers));
}

TEST(KmerSetCompact, Size) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 100000;

  const KmerSet<K, N, KeyType> kmer_set =
      GetRandomKmerSet<K, N, KeyType>(n, true);

  const KmerSetCompact<K, N, KeyType> kmer_set_compact =
      KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, true, true,
                                                 n_workers);

  ASSERT_EQ(kmer_set.Size(), kmer_set_compact.Size());
}

TEST(KmerSetCompact, FromAndToKmerSet) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 100000;

  const KmerSet<K, N, KeyType> kmer_set =
      GetRandomKmerSet<K, N, KeyType>(n, true);

  const KmerSetCompact<K, N, KeyType> kmer_set_compact =
      KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, true, true,
                                                 n_workers);

  const KmerSet<K, N, KeyType> reconstructed =
      kmer_set_compact.ToKmerSet(true, n_workers);

  ASSERT_TRUE(kmer_set.Equals(reconstructed, n_workers));
}

TEST(KmerSetCompact, GetBloomFilter) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 100000;

  const KmerSet<K, N, KeyType> kmer_set =
      GetRandomKmerSet<K, N, KeyType>(n, true);

  const KmerSetCompact<K, N, KeyType> kmer_set_compact =
      KmerSetCompact<K, N, KeyType>::FromKmerSet(kmer_set, true, true,
                                                 n_workers);

  const std::vector<bool> bloom_filter =
      kmer_set_compact.GetBloomFilter(n, n * 10, n_workers);

  ASSERT_EQ(bloom_filter.size(), n);
  ASSERT_FALSE(std::vector<bool>(n) == bloom_filter);
}
