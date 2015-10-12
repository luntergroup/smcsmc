Record recombination opportunity before extend to the new genealogy

##Package status

package  | Build Status
-------- | -----------------
smcsmc   | [![Build Status](https://magnum.travis-ci.com/luntergroup/smcsmc.svg?branch=master)](https://magnum.travis-ci.com/luntergroup/smcsmc)
scrm     | [![Build Status](https://travis-ci.org/luntergroup/scrm_jz_stable_branch.svg?branch=smcsmcSCRM)](https://travis-ci.org/luntergroup/scrm_jz_stable_branch)


# RECOMBINATION CALCULATION IS CURRENTLY OFF ...
build-dep libcppunit-dev


## Installation
On Debian/Ubuntu based systems:
```bash
apt-get install build-essential autoconf autoconf-archive doxygen libcppunit-dev
CXXFLAGS="-O3" ./bootstrap
make
```

On Mac OS:
```bash
port install autoconf autoconf-archive doxygen cppunit
./bootstrap
make
```

Still working progress
In order to use the multithread version on Mac, please follow the instruction
build OpenMp on Mac (see http://clang-omp.github.io)
first download
$ git clone https://github.com/clang-omp/llvm
$ git clone https://github.com/clang-omp/compiler-rt llvm/projects/compiler-rt
$ git clone -b clang-omp https://github.com/clang-omp/clang llvm/tools/clang




29-Jun-14

clang is used by mac by default, if want to use multithreading version, should suggest to use the gcc compiler.

no longer requires boost ...

