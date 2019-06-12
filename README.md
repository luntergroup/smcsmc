## Demographic Inference using Particle Filters

SMCSMC (Sequential Monte Carlo for the Sequential Markovian Coalescent) or SMC2 is a program for inferring population history from multiple genome sequences. 

## Installation

This repository contains two effectively independant components, and both must be installed to properly use SMC2. We have automated this process in a `conda` package, and we highly recommend installing it this way.

```sh
conda install -c terhorst, conda-forge smcsmc
```

## Running

```./smc2
```

For examples see the documentation in `docs/build/html/` or [online](https://github.com/luntergroup/smcsmc/tree/master/docs/build/html/)

### Installation from Source

Alternatively, a combination of `cmake` and `pip` can be used to install the python and core components.

#### To install dependencies

Download, or use a package manager, to install the following packages:

- boost
- cmake
- tcmalloc

#### To install the backend

```sh
mkdir build; cd build
cmake ..
cmake
```

#### To install the frontend

```sh
pip install -r dependencies
pip install .
```
