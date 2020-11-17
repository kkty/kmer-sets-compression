#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "kmer.h"
#include "kmer_set.h"
#include "kmer_set_compressed.h"
#include "spdlog/spdlog.h"

ABSL_FLAG(bool, debug, false, "enable debugging messages");
ABSL_FLAG(std::string, decompressor, "cat", "specify decompressor");
ABSL_FLAG(int, workers, 1, "number of workers");
ABSL_FLAG(bool, canonical, false, "use canonical k-mers");

int main(int argc, char** argv) {
  const std::string file_name = absl::ParseCommandLine(argc, argv)[1];

  if (absl::GetFlag(FLAGS_debug)) spdlog::set_level(spdlog::level::debug);

  const int K = 17;
  const int n_workers = absl::GetFlag(FLAGS_workers);
  using KeyType = uint16_t;

  spdlog::info("constructing kmer_set_compressed");
  const KmerSetCompressed<K, KeyType> kmer_set_compressed =
      KmerSetCompressed<K, KeyType>::Load(
          file_name, absl::GetFlag(FLAGS_decompressor));
  spdlog::info("constructed kmer_set_compressed");

  spdlog::info("constructing kmer_set");
  const KmerSet<K, KeyType> kmer_set =
      kmer_set_compressed.ToKmerSet(absl::GetFlag(FLAGS_canonical), n_workers);
  spdlog::info("constructed kmer_set");

  spdlog::info("kmer_set.Size() = {}", kmer_set.Size());
}
