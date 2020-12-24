## Requirements

- Unix-like system
- C++ compiler with C++17 support
- [CMake](https://cmake.org/)

## Installing dependencies

The following script installs [Abseil](https://abseil.io/), [Boost](https://www.boost.org/)
, [Googletest](https://github.com/google/googletest), [mimalloc](https://github.com/microsoft/mimalloc)
, [spdlog](https://github.com/gabime/spdlog), and [streamvbyte](https://github.com/lemire/streamvbyte) in `extern`
directory.

```shell
cd extern
PREFIX=$(pwd) ./install.sh
```

## Building

The following script builds the programs in `build` direcory.

```shell
mkdir build
cd build
cmake .. -DCMAKE_PREFIX_PATH=$(pwd)/../extern -DCMAKE_BUILD_TYPE=Release
make -j
```

## Programs

```
$ ./build/kmerset-compress
kmerset-compress: Reads a FASTA file and constructs a set of k-mers. Usage: ./build/kmerset-compress [options] <path to file>

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --check (does compression & decompression to see if it is working
      correctly); default: false;
    --compressor (an program to compress output files; e.g., "bzip2" for bzip2,
      "" for no compression); default: "";
    --cutoff (ignore k-mers that appear less often than this value); default: 1;
    --debug (enable debugging messages); default: false;
    --decompressor (an program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "" for no decompression); default: "";
    --k (the length of k-mers); default: 15;
    --out (output file name); default: "";
    --workers (number of threads to use); default: 1;
```

```
$ ./build/kmerset-decompress
kmerset-decompress: Prints the metadata of a k-mer set. Usage: ./build/kmerset-decompress [options] <path to file>

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --debug (enable debugging messages); default: false;
    --decompressor (an program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "" for no decompression); default: "";
    --k (the length of k-mers); default: 15;
    --workers (number of threads to use); default: 1;
```

```
$ ./build/kmersetset-compress --help
kmersetset-compress: Compresses multiple k-mer sets. Usage: ./build/kmersetset-compress [options] <paths to file> <path to file> ...

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --compressor (program to compress dumped file); default: "";
    --debug (enable debugging messages); default: false;
    --decompressor (an program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "" for no decompression); default: "";
    --extension (extension for output files); default: "bin";
    --iteration (number of iterations for KmerSetSet); default: 1;
    --k (the length of k-mers); default: 15;
    --out (directory path to save dumped files); default: "";
    --out_graph (path to save dumped DOT file); default: "";
    --parallel_input (read files in parallel); default: false;
    --workers (number of threads to use); default: 1;
```

```
$ ./build/kmersetset-decompress --help
kmersetset-decompress: Decompresses the output of "kmerset-compress". Usage: ./build/kmersetset-decompress [options] <path to directory>

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --debug (enable debugging messages); default: false;
    --decompressor (an program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "" for no decompression); default: "";
    --extension (extension of files in folder); default: "txt";
    --k (the length of k-mers); default: 15;
    --workers (number of threads to use); default: 1;
```

```
$ ./build/kmerset-stat
kmerset-stat: For each pair in multiple k-mer sets, calculates their similarity. Usage: ./build/kmerset-stat [options] <paths to file> <path to file> ...

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --debug (enable debugging messages); default: false;
    --decompressor (an program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "" for no decompression); default: "";
    --k (the length of k-mers); default: 15;
    --workers (number of threads to use); default: 1;
```

```
$ ./build/kmerset-immutable
kmerset-immutable: Runs a benchmark for intersecting 2 k-mer sets. Usage: ./build/kmerset-immutable [options] <path to file> <path to file>

  Flags:
    --canonical (set this flag when handling canonical k-mers); default: true;
    --debug (enable debugging messages); default: false;
    --decompressor (an program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "" for no decompression); default: "";
    --k (the length of k-mers); default: 15;
    --repeats (number of repeats); default: 100;
    --workers (number of threads to use); default: 1;
```

```
$ ./build/spss --help
spss: Runs a benchmark for SPSS construction using a single k-mer set. Usage: ./build/spss [options] <path to file>

  Flags:
    --buckets (number of buckets for SPSS calculation); default: 1;
    --debug (enable debugging messages); default: false;
    --decompressor (an program to decompress input files; e.g., "bzip2 -d" for
      bzip2, "" for no decompression); default: "";
    --k (the length of k-mers); default: 15;
    --repeats (number of repeats); default: 1;
    --workers (number of threads to use); default: 1;
```