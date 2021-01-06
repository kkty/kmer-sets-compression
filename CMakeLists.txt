cmake_minimum_required(VERSION 3.10)

project(research)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_LINKER_FLAGS_DEBUG "-fsanitize=address")

find_package(absl REQUIRED)
find_package(GTest REQUIRED)
find_package(spdlog REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread)
find_package(mimalloc REQUIRED)
find_package(benchmark REQUIRED)

find_library(streamvbyte_LIBRARIES streamvbyte_static REQUIRED)

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
        absl::statusor
        absl::random_random
        ${Boost_LIBRARIES}
        ${streamvbyte_LIBRARIES}
)

foreach (
        NAME
        kmerset-build
        kmerset-intersection-approximation
        kmerset-intersection-performance
        kmerset-multiple-compress
        kmerset-multiple-decompress
        kmerset-similarity
        kmerset-stat
        spss
)
    add_executable(${NAME} src/${NAME}.cc)
endforeach ()

include(GoogleTest)
include(CTest)

foreach (
        NAME
        int_set
        kmer
        kmer_counter
        kmer_set
        kmer_set_immutable
        kmer_set_set
        parallel_disjoint_set
        range
        io
        spss
)
    add_executable(test_${NAME} test/${NAME}.cc)
    target_link_libraries(test_${NAME} GTest::gtest_main)
    gtest_add_tests(test_${NAME} "" AUTO)
endforeach ()

foreach (
        NAME
        int_set
        kmer
        kmer_counter
        kmer_set
        spss
)
    add_executable(benchmark_${NAME} benchmark/${NAME}.cc)
    target_link_libraries(benchmark_${NAME} benchmark::benchmark)
endforeach()
