cmake_minimum_required (VERSION 3.0)

project(research)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_LINKER_FLAGS_DEBUG "-fno-omit-frame-pointer -fsanitize=address")

find_package(absl REQUIRED)
find_package(GTest REQUIRED)
find_package(spdlog REQUIRED)
find_package(Boost REQUIRED)
find_package(mimalloc REQUIRED)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
include_directories(lib)

link_libraries(
  mimalloc
  spdlog::spdlog
  absl::strings
  absl::flat_hash_set
  absl::flat_hash_map
  absl::flags
  absl::flags_parse
  absl::status
)

add_executable(search src/search.cc)
add_executable(compact src/compact.cc)
add_executable(compress src/compress.cc)

link_libraries(
  GTest::gtest
  GTest::gtest_main
)

add_executable(test_kmer test/kmer.cc)
add_executable(test_kmer_set test/kmer_set.cc)
add_executable(test_kmer_set_compact test/kmer_set_compact.cc)
add_executable(test_suffix_array test/suffix_array.cc)
add_executable(test_range test/range.cc)
add_executable(test_graph test/graph.cc)
