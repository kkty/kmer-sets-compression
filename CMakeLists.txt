cmake_minimum_required(VERSION 3.10)

project(research)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_LINKER_FLAGS_DEBUG "-fsanitize=address")

find_package(absl REQUIRED)
find_package(GTest REQUIRED)
find_package(spdlog REQUIRED)
find_package(Boost REQUIRED COMPONENTS system thread)
find_package(mimalloc REQUIRED)

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
        kmerset-compress
        kmerset-decompress
        kmerset-stat
        kmersetset-compress
        kmersetset-decompress
        search
        search-allpair
        spss
)
    add_executable(${NAME} src/${NAME}.cc)
endforeach ()

link_libraries(GTest::gtest_main)
include(GoogleTest)
include(CTest)

foreach (
        NAME
        int_set
        kmer
        kmer_counter
        kmer_set_compact
        kmer_set
        range
        io
        search
        unitigs
)
    add_executable(test_${NAME} test/${NAME}.cc)
    gtest_add_tests(test_${NAME} "" AUTO)
endforeach ()
