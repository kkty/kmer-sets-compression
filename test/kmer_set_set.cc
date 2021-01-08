#include "core/kmer_set_set.h"

#include <cstdint>
#include <string>
#include <utility>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "core/kmer_set.h"
#include "core/kmer_set_compact.h"
#include "gtest/gtest.h"
#include "io.h"
#include "random.h"

TEST(KmerSetSet, GetRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 10;
  const int m = 10000;

  std::vector<KmerSetCompact<K, N, KeyType>> kmer_sets_compact =
      GetRandomKmerSetsCompact<K, N, KeyType>(n, m, true, n_workers);

  const KmerSetSet<K, N, KeyType> kmer_set_set(kmer_sets_compact, true,
                                               n_workers);

  for (int i = 0; i < n; i++) {
    ASSERT_TRUE(kmer_set_set.Get(i, true, n_workers)
                    .Equals(kmer_sets_compact[i].ToKmerSet(true, n_workers),
                            n_workers));
  }
}

TEST(KmerSetSet, DumpAndLoadRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 10;
  const int m = 10000;

  const TemporaryDirectory temporary_directory;

  KmerSetSet<K, N, KeyType> kmer_set_set =
      GetRandomKmerSetSet<K, N, KeyType>(n, m, true, n_workers);

  {
    absl::Status status =
        kmer_set_set.Dump(temporary_directory.Name(), "", "txt", n_workers);
    ASSERT_TRUE(status.ok());
  }

  KmerSetSet<K, N, KeyType> loaded;
  {
    absl::StatusOr<KmerSetSet<K, N, KeyType>> statusor =
        KmerSetSet<K, N, KeyType>::Load(temporary_directory.Name(), "", "txt",
                                        n_workers);

    ASSERT_TRUE(statusor.ok());

    loaded = std::move(statusor).value();
  }

  ASSERT_EQ(kmer_set_set.Size(), loaded.Size());

  // Makes sure that every kmer set can be reconstructed.

  for (int i = 0; i < kmer_set_set.Size(); i++) {
    ASSERT_TRUE(kmer_set_set.Get(i, true, n_workers)
                    .Equals(loaded.Get(i, true, n_workers), n_workers));
  }
}

TEST(KmerSetSet, ReaderRandom) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 10;
  const int m = 10000;

  const TemporaryDirectory temporary_directory;

  KmerSetSet<K, N, KeyType> kmer_set_set =
      GetRandomKmerSetSet<K, N, KeyType>(n, m, true, n_workers);

  {
    absl::Status status =
        kmer_set_set.Dump(temporary_directory.Name(), "", "txt", n_workers);
    ASSERT_TRUE(status.ok());
  }

  KmerSetSetReader<K, N, KeyType> kmer_set_set_reader;
  {
    absl::StatusOr<KmerSetSetReader<K, N, KeyType>> statusor =
        KmerSetSetReader<K, N, KeyType>::FromDirectory(
            temporary_directory.Name(), "txt", "", true);

    ASSERT_TRUE(statusor.ok());

    kmer_set_set_reader = std::move(statusor).value();
  }

  ASSERT_EQ(kmer_set_set.Size(), kmer_set_set_reader.Size());

  // Makes sure that every kmer set can be reconstructed.

  for (int i = 0; i < kmer_set_set.Size(); i++) {
    absl::StatusOr<KmerSet<K, N, KeyType>> statusor =
        kmer_set_set_reader.Get(i, n_workers);

    ASSERT_TRUE(statusor.ok());

    const KmerSet<K, N, KeyType> kmer_set = std::move(statusor).value();

    ASSERT_TRUE(kmer_set_set.Get(i, true, n_workers).Equals(kmer_set, n_workers));
  }
}
