#!/bin/bash

# PREFIX should be set as an environment variable.

set -eu

WORK_DIR=$(mktemp -d)

echo "working directory: $WORK_DIR"

pushd $WORK_DIR
  echo "installing boost"
  wget https://dl.bintray.com/boostorg/release/1.75.0/source/boost_1_75_0.tar.gz
  tar xf boost_1_75_0.tar.gz
  pushd boost_1_75_0
    ./bootstrap.sh --prefix=$PREFIX --with-libraries=thread
    ./b2 install -j4
  popd

  echo "installing abseil"
  wget https://github.com/abseil/abseil-cpp/archive/20200923.2.tar.gz
  tar xf 20200923.2.tar.gz
  pushd abseil-cpp-20200923.2
    mkdir -p build
    pushd build
      cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=$PREFIX
      make -j
      make install
    popd
  popd

  echo "installing googletest"
  wget https://github.com/google/googletest/archive/release-1.10.0.tar.gz
  tar xf release-1.10.0.tar.gz
  pushd googletest-release-1.10.0
    mkdir -p build
    pushd build
      cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX
      make -j
      make install
    popd
  popd

  echo "installing mimalloc"
  wget https://github.com/microsoft/mimalloc/archive/v1.6.7.tar.gz
  tar xf v1.6.7.tar.gz
  pushd mimalloc-1.6.7
    mkdir build
    pushd build
      cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX
      make -j
      make install
    popd
  popd

  echo "installing spdlog"
  wget https://github.com/gabime/spdlog/archive/v1.8.2.tar.gz
  tar xf v1.8.2.tar.gz
  pushd spdlog-1.8.2
    mkdir build
    pushd build
      cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX
      make -j
      make install
    popd
  popd

  echo "installing google benchmark"
  wget https://github.com/google/benchmark/archive/v1.5.2.tar.gz
  tar xf v1.5.2.tar.gz
  pushd benchmark-1.5.2
    mkdir build
    pushd build
      cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DBENCHMARK_DOWNLOAD_DEPENDENCIES=ON
      make -j
      make install
    popd
  popd
popd

rm -rf $WORK_DIR