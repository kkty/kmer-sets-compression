cmake_minimum_required (VERSION 3.0)

project(search)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fvisibility=hidden")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_LINKER_FLAGS_DEBUG "-fno-omit-frame-pointer -fsanitize=address")

add_subdirectory(extern/abseil-cpp)
add_subdirectory(extern/mimalloc)
add_subdirectory(extern/spdlog)
add_subdirectory(extern/googletest)

find_package(Boost REQUIRED iostreams)

add_executable(search src/search.cc)
add_executable(compact src/compact.cc)
add_executable(compress src/compress.cc)

add_executable(test_kmer test/kmer.cc)
add_executable(test_kmer_set test/kmer_set.cc)
add_executable(test_kmer_set_compact test/kmer_set_compact.cc)
add_executable(test_suffix_array test/suffix_array.cc)
add_executable(test_range test/range.cc)
add_executable(test_graph test/graph.cc)

target_include_directories(search PUBLIC lib)
target_include_directories(compact PUBLIC lib)
target_include_directories(compress PUBLIC lib)

target_include_directories(test_kmer PUBLIC lib)
target_include_directories(test_kmer_set PUBLIC lib)
target_include_directories(test_suffix_array PUBLIC lib)
target_include_directories(test_kmer_set_compact PUBLIC lib)
target_include_directories(test_range PUBLIC lib)
target_include_directories(test_graph PUBLIC lib)

set(
  src_libs
  mimalloc
  spdlog
  absl::strings
  absl::flat_hash_set
  absl::flat_hash_map
  absl::flags
  absl::flags_parse
  ${Boost_LIBRARIES}
)

set (
  test_libs
  gtest
  gtest_main
  ${src_libs}
)

target_link_libraries(search ${src_libs})
target_link_libraries(compact ${src_libs})
target_link_libraries(compress ${src_libs})

target_link_libraries(test_kmer ${test_libs})
target_link_libraries(test_kmer_set ${test_libs})
target_link_libraries(test_suffix_array ${test_libs})
target_link_libraries(test_kmer_set_compact ${test_libs})
target_link_libraries(test_range ${test_libs})
target_link_libraries(test_graph ${test_libs})
