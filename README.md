# Demographic Inference using a Particle Filter
[![Anaconda-Server Badge](https://anaconda.org/anaconda/anaconda/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)

SMCSMC (Sequential Monte Carlo for the Sequential Markovian Coalescent) or SMC2 is a program for inferring population history from multiple genome sequences. It includes both a python package `smcsmc` and a command line interface `smc2` along with two backend binaries `smcsmc`/`scrm`.

For examples and explaination, please see the documentation in `docs/`. 

## Installation

This repository contains two components, and both must be installed to properly use `smcsmc`.

### Recommended Installation via `conda`

We have automated this process in a `conda` package, and we highly recommend installing it this way.

```sh
conda install -c terhorst -c conda-forge smcsmc
```

> **We are in the process of uploading our package to `conda` and will remove this message when it is available. Until then, please use the manual installation method given below.**

### Installation from Source

Alternatively, a combination of `cmake` and `pip` can be used to install the python and core components.

#### Obtain the code

```sh
git clone git@github.com:luntergroup/smcsmc.git git-smcsmc
cd git-smcsmc
git submodule init
git submodule update
```

#### Install dependencies

Download, or use a package manager, to install the following packages:

- boost
- cmake
- tcmalloc

#### Install the c++ backend

```sh
mkdir build; cd build
cmake ..
cmake
```

#### Install the frontend

```sh
pip install -r dependencies
pip install .
```

## Citation

If you use `smcsmc` in your work, please cite the following articles:

1. Henderson, D., Zhu, S. (Joe), & Lunter, G. (2018). Demographic inference using particle filters for continuous Markov jump processes. BioRxiv, 382218. https://doi.org/10.1101/382218

2. Staab, P. R., Zhu, S., Metzler, D., & Lunter, G. (2015). scrm: efficiently simulating long sequences using the approximated coalescent with recombination. Bioinformatics, 31(10), 1680–1682. https://doi.org/10.1093/bioinformatics/btu861


