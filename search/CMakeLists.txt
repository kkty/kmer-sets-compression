cmake_minimum_required (VERSION 3.0)

project(search)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(extern/abseil-cpp)
add_subdirectory(extern/mimalloc)
add_subdirectory(extern/spdlog)
add_subdirectory(extern/googletest)

find_package(OpenMP REQUIRED)

add_executable(main main.cc)
add_executable(test test.cc)

target_include_directories(main PUBLIC lib)
target_include_directories(test PUBLIC lib)

target_link_libraries(
    test
    mimalloc
    spdlog
    absl::strings
    absl::flat_hash_set
    absl::flat_hash_map
    absl::flags
    absl::flags_parse
    gtest
    gtest_main
)

target_link_libraries(
    main
    mimalloc
    spdlog
    absl::strings
    absl::flat_hash_set
    absl::flat_hash_map
    absl::flags
    absl::flags_parse
)
