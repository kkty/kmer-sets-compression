```
# On Ubuntu20.04

sudo apt update
sudo apt install gcc-10 g++-10 gcc g++

# gperftools

sudo apt install autoconf libtool make
git clone https://github.com/gperftools/gperftools
pushd gperftools
./autogen.sh
./configure
make -j
sudo make install
popd

# sparsehash

git clone https://github.com/sparsehash/sparsehash
pushd sparsehash
./configure
make -j
sudo make install
popd

# add "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib" in .bashrc

g++-10 -std=c++17 -ltcmalloc -O3 -fopenmp main-compact.cc

./a.out SRR957915.fastq
```
