cmake_minimum_required (VERSION 3.0)

project(20201024)


set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(extern/abseil-cpp)
add_subdirectory(extern/mimalloc)
add_subdirectory(extern/spdlog)

find_package(OpenMP REQUIRED)

add_executable(main main.cc)
target_link_libraries(main absl::strings absl::flat_hash_set mimalloc OpenMP::OpenMP_CXX spdlog)

