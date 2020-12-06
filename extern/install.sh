#!/bin/bash

# PREFIX should be set as an environment variable.

set -eu

for DIR in abseil-cpp googletest mimalloc spdlog streamvbyte
do
  pushd $DIR
    mkdir -p build
    pushd build
      cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_FLAGS_RELEASE='-O3 -march=native' -DCMAKE_INSTALL_PREFIX=$PREFIX
      make -j
      make install
    popd
    rm -rf build
  popd
done