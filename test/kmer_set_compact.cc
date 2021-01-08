#include "core/kmer_set_compact.h"

#include <algorithm>
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

TEST(KmerSetCompact, GetSampledKmerSet) {
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

  std::vector<int> bucket_ids;
  for (int i = 0; i < (1 << N); i++) bucket_ids.push_back(i);

  std::reverse(bucket_ids.begin(), bucket_ids.end());

  KmerSet<K, N, KeyType> reconstructed;

  {
    std::vector<std::vector<KeyType>> sampled =
        kmer_set_compact.GetSampledKmerSet(bucket_ids, true, n_workers);

    for (int i = 0; i < (1 << N); i++) {
      ASSERT_TRUE(std::is_sorted(sampled[i].begin(), sampled[i].end()));

      for (KeyType key : sampled[i]) {
        reconstructed.Add(
            GetKmerFromBucketAndKey<K, N, KeyType>(bucket_ids[i], key));
      }
    }
  }

  ASSERT_TRUE(kmer_set.Equals(reconstructed, n_workers));
}