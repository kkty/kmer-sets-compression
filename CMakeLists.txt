cmake_minimum_required (VERSION 3.0)

project(research)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_LINKER_FLAGS_DEBUG "-fsanitize=address")

include(CTest)

find_package(absl REQUIRED)
find_package(GTest REQUIRED)
find_package(spdlog REQUIRED)
find_package(Boost REQUIRED COMPONENTS system thread)
find_package(mimalloc REQUIRED)
find_package(LEMON REQUIRED)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${LEMON_INCLUDE_DIRS})
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
  absl::statusor
  absl::random_random
  ${Boost_LIBRARIES}
  ${LEMON_LIBRARIES}
)

foreach(
  NAME
  search
  compress
  decompress
  kmersetset
)
  add_executable(${NAME} src/${NAME}.cc)
endforeach()

link_libraries(GTest::gtest_main)

foreach(
  NAME
  kmer
  kmer_counter
  kmer_set
  range
  graph
  neighbor_joining
  io
  search
)
  add_executable(test_${NAME} test/${NAME}.cc)
  add_test(test_${NAME} test_${NAME})
endforeach()

