#include "core/kmer_set_set.h"

#include <cstdint>
#include <filesystem>
#include <ios>
#include <limits>
#include <sstream>
#include <string>

#include "absl/random/random.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "core/kmer_set.h"
#include "core/kmer_set_immutable.h"
#include "gtest/gtest.h"
#include "random.h"

class TemporaryDirectory {
 public:
  TemporaryDirectory() {
    std::stringstream ss;
    absl::InsecureBitGen bitgen;

    ss << std::hex
       << absl::Uniform<std::uint64_t>(
              absl::IntervalClosed, bitgen, 0,
              std::numeric_limits<std::uint64_t>::max());

    name_ = (std::filesystem::temp_directory_path() / ss.str()).string();
  }

  ~TemporaryDirectory() { std::filesystem::remove_all(name_); }

  std::string Name() const { return name_; }

 private:
  std::string name_;
};

TEST(KmerSetSet, CtorAndGet) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 10;
  const int m = 10000;

  std::vector<KmerSetImmutable<K, N, KeyType>> kmer_sets_immutable =
      GetRandomKmerSetsImmutable<K, N, KeyType>(n, m, n_workers);

  const KmerSetSet<K, N, KeyType> kmer_set_set(kmer_sets_immutable, true,
                                               n_workers);

  for (int i = 0; i < n; i++) {
    ASSERT_TRUE(
        kmer_set_set.Get(i, n_workers)
            .Equals(kmer_sets_immutable[i].ToKmerSet(n_workers), n_workers));
  }
}

TEST(KmerSetSet, DumpAndLoad) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 10;
  const int m = 10000;

  const TemporaryDirectory temporary_directory;

  KmerSetSet<K, N, KeyType> kmer_set_set =
      GetRandomKmerSetSet<K, N, KeyType>(n, m, n_workers);

  {
    absl::Status status = kmer_set_set.Dump(temporary_directory.Name(), "",
                                            "txt", true, false, n_workers);
    ASSERT_TRUE(status.ok());
  }

  KmerSetSet<K, N, KeyType> loaded;
  {
    absl::StatusOr<KmerSetSet<K, N, KeyType>> statusor =
        KmerSetSet<K, N, KeyType>::Load(temporary_directory.Name(), "", "txt",
                                        true, n_workers);

    ASSERT_TRUE(statusor.ok());

    loaded = std::move(statusor).value();
  }

  ASSERT_EQ(kmer_set_set.Size(), loaded.Size());

  // Makes sure that every kmer set can be reconstructed.

  for (int i = 0; i < kmer_set_set.Size(); i++) {
    ASSERT_TRUE(kmer_set_set.Get(i, n_workers)
                    .Equals(loaded.Get(i, n_workers), n_workers));
  }
}

TEST(KmerSetSet, Reader) {
  const int K = 9;
  const int N = 10;
  using KeyType = std::uint8_t;
  const int n_workers = 4;

  const int n = 10;
  const int m = 10000;

  const TemporaryDirectory temporary_directory;

  KmerSetSet<K, N, KeyType> kmer_set_set =
      GetRandomKmerSetSet<K, N, KeyType>(n, m, n_workers);

  {
    absl::Status status = kmer_set_set.Dump(temporary_directory.Name(), "",
                                            "txt", true, false, n_workers);
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

    ASSERT_TRUE(kmer_set_set.Get(i, n_workers).Equals(kmer_set, n_workers));
  }
}
