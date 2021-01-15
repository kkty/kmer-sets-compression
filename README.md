# Compression of Multiple k-mer Sets

## Requirements

- Unix-like system
- C++ compiler with C++17 support
- [CMake](https://cmake.org/) (>= 3.17)

### Tested Configurations

- Ubuntu 20.04, GCC 10, and CMake 3.19
- MacOS Big Sur, Clang 12, and CMake 3.18

## Installing Dependencies

The following script installs [Abseil](https://abseil.io/),
[Boost](https://www.boost.org/),
[GoogleTest](https://github.com/google/googletest),
[Google Benchmark](https://github.com/google/benchmark),
[mimalloc](https://github.com/microsoft/mimalloc),
[streamvbyte](https://github.com/lemire/streamvbyte), and [spdlog](https://github.com/gabime/spdlog) in `extern`
directory.

```shell
cd extern
PREFIX=$(pwd) ./install.sh
```

## Building

The following script builds programs in `build` direcory.

```shell
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=$(pwd)/../extern -DCMAKE_BUILD_TYPE=Release
make -j
```

## Programs

The following executables are built in the build directory.

### kmerset-build

```
$ ./kmerset-build --help     
kmerset-build: Reads a FASTA file and constructs a set of k-mers. Usage: ./kmerset-build [options] <path to file>

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --check (does compression & decompression to see if it is working
      correctly); default: false;
    --compressor (a program to compress output files; e.g., "bzip2" for bzip2,
      "gzip" for gzip, and "" for no compression); default: "";
    --cutoff (ignore k-mers that appear less often than this value); default: 1;
    --debug (enable debugging messages); default: false;
    --decompressor (a program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "gzip -d" for gzip, and "" for no decompression); default: "";
    --k (the length of k-mers); default: 15;
    --out (output file name); default: "";
    --workers (number of threads to use); default: 1;

Try --helpfull to get a list of all flags.
```

#### Example

The following command reads `foo.fasta.gz`, counts canonical k-mers, removes ones that appear less than 4 times, and
saves the resulting k-mer set data to `foo.kmerset.bz2`. k is set to 23. 8 threads will be used.

```
./kmerset-build --canonical --compressor='bzip2' --cutoff=4 --decompressor='gzip2 -d' --k=23 --out=foo.kmerset.bz2 --workers=8 foo.fasta.gz
```

### kmerset-stat

```
$ ./kmerset-stat --help 
kmerset-stat: Prints the metadata of a k-mer set. Usage: ./kmerset-stat [options] <path to file>

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --debug (enable debugging messages); default: false;
    --decompressor (a program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "gzip -d" for gzip, and "" for no decompression); default: "";
    --k (the length of k-mers); default: 15;
    --workers (number of threads to use); default: 1;

Try --helpfull to get a list of all flags.
```

#### Example

The following command shows the metadata of the k-mer set represented by `foo.kmerset.bz2`. It is assumed that the file
represents a k-mer set of canonical k-mers where k is 23. 8 threads will be used.

```
./kmerset-stat --canonical --decompressor='bzip2 -d' --k=23 --workers=8 foo.kmerset.bz2
```

### kmerset-multiple-compress

```
$ ./kmerset-multiple-compress --help 
kmerset-multiple-compress: Compresses multiple k-mer sets. Usage: ./kmerset-multiple-compress [options] <paths to file> <path to file> ...

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --compressor (a program to compress output files; e.g., "bzip2" for bzip2,
      "gzip" for gzip, and "" for no compression); default: "";
    --debug (enable debugging messages); default: false;
    --decompressor (a program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "gzip -d" for gzip, and "" for no decompression); default: "";
    --extension (extension for output files); default: "txt";
    --k (the length of k-mers); default: 15;
    --out (directory path to save dumped files); default: "";
    --out_graph (path to save dumped DOT file); default: "";
    --parallel_input (read files in parallel); default: false;
    --workers (number of threads to use); default: 1;

Try --helpfull to get a list of all flags.
```

#### Example

The following command reads `./data/*.kmerset.gz`, and compresses the obtained k-mer sets. The output will be saved
to `./compressed/*.bz2` after bzip2-ed. The DOT file representing the graph data will be saved to `./graph.gv`. It
handles canonical k-mers where k is 23. 8 threads will be used.

```
./kmer-set-multiple-compress --canonical --compressor='bzip2' --decompressor='gzip -d' --extension='bz2' --k=23 --out=./compressed --out_graph=./graph.gv --workers=8 ./data/*.kmerset.gz
```

### kmerset-multiple-decompress

```
$ ./kmerset-multiple-decompress --help
kmerset-multiple-decompress: Decompresses the output of "kmerset-multiple-compress". Usage: ./kmerset-multiple-decompress [options] <path to directory>

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --debug (enable debugging messages); default: false;
    --decompressor (a program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "gzip -d" for gzip, and "" for no decompression); default: "";
    --extension (extension of files in folder); default: "txt";
    --k (the length of k-mers); default: 15;
    --workers (number of threads to use); default: 1;

Try --helpfull to get a list of all flags.
```

#### Example

The following command reads the output of `kmerset-multiple-compress` saved to `./compressed/*.bz2`, decompresses the
compressed data, and prints the metadata of each of the original k-mer sets. It handles canonical k-mers where k is 23.
8 threads will be used.

```
./kmerset-multiple-decompress --canonical --decompressor='bzip2 -d' --extension='bz2' --k=23 --workers=8 ./compressed
```

### spss-benchmark

```
$ ./spss-benchmark --help
spss-benchmark: Runs a benchmark for SPSS construction using a single k-mer set. Usage: ./spss-benchmark [options] <path to file>

  Flags from Users/kazushi/work/research/src/spss-benchmark.cc:
    --buckets (number of buckets for SPSS calculation); default: 1;
    --debug (enable debugging messages); default: false;
    --decompressor (a program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "gzip -d" for gzip, and "" for no decompression); default: "";
    --k (the length of k-mers); default: 15;
    --repeats (number of repeats); default: 1;
    --workers (number of threads to use); default: 1;

Try --helpfull to get a list of all flags.
```

#### Example

The following command loads the k-mer set represented by `foo.kmerset.bz2`, and runs a benchmark to compare our proposed
SPSS construction algorithm with UST algorithm. It will use the k value of 23. 1024 buckets (a parameter for the
propsoed algorithm) and 8 threads will be used.

```
./spss-benchmark --buckets=1024 --decompressor='bzip2 -d' --k=23 --repeats=10 --workers=8 foo.kmerset.bz2
```

#### Note

The input file should contain canonical k-mers.

## Running Tests

The following command, when executed in the build directory, invokes all the tests.

```
ctest
```

It is also possible to configure test execution by providing arguments to `ctest`. Refer
to [ctest documentation](https://cmake.org/cmake/help/latest/manual/ctest.1.html) for details.

## Project Structure

- `lib` contains most of the source code. The code in `lib/core` provides core functionalities, and the code outside the
  directory provides helper functions. The files in `lib/core` do not depend on the files outside the `lib/core`
  directory.
- `src` contains source codes for executables. Each `.cc` file corresponds to one executable with the same name.
- `test` contains source code for functions and classes defined in `lib/core`.
- `benchmark` contains source code for benchmarks for critical functions and classes.