# A Particle Filter for Demographic Inference
[![Anaconda-Server Badge](https://anaconda.org/luntergroup/smcsmc/badges/version.svg)](https://anaconda.org/luntergroup/smcsmc) [![Anaconda-Server Badge](https://anaconda.org/luntergroup/smcsmc/badges/platforms.svg)](https://anaconda.org/luntergroup/smcsmc) [![Documentation Status](https://readthedocs.org/projects/smcsmc/badge/?version=latest)](https://smcsmc.readthedocs.io/en/latest/?badge=latest) [![Anaconda-Server Badge](https://anaconda.org/luntergroup/smcsmc/badges/downloads.svg)](https://anaconda.org/luntergroup/smcsmc)
 

SMCSMC (Sequential Monte Carlo for the Sequential Markovian Coalescent) or SMC2 is a program for inferring population history from multiple genome sequences. It includes both a python package `smcsmc` and a command line interface `smc2` along with two backend binaries `smcsmc`/`scrm`.

For examples and explaination, please see the documentation [online](https://smcsmc.readthedocs.io) or in `docs/`.

## Installation

This repository contains two components, and both must be installed to properly use `smcsmc`.

### Recommended Installation via `conda`

We have automated this process in a `conda` package, and we highly recommend installing it this way.

```sh
conda config --add channels conda-forge
conda config --add channels terhorst
conda install -c luntergroup smcsmc
```

> NOTE: We currently only support `conda` installation on **64 bit Linux**.  If you are using a different operating system you must install manually -- see below. 

### Manual installation

Alternatively, a combination of `cmake` and `pip` can be used to install the python and core components:

#### Obtain the code

```sh
git clone git@github.com:luntergroup/smcsmc.git git-smcsmc
cd git-smcsmc
git submodule init
git submodule update
```

#### Install dependencies

Download and install the following packages (or use a package manager):

- boost
- cmake
- tcmalloc

#### Install the c++ backend

```sh
mkdir build; cd build
cmake ..
make
cd ..
```

#### Install the frontend

```sh
pip install -r requirements.txt
pip install .
```

## Citation

If you use `smcsmc` in your work, please cite the following articles:

1. Henderson, D., Zhu, S. (Joe), & Lunter, G. (2018). Demographic inference using particle filters for continuous Markov jump processes. BioRxiv, 382218. https://doi.org/10.1101/382218

2. Staab, P. R., Zhu, S., Metzler, D., & Lunter, G. (2015). scrm: efficiently simulating long sequences using the approximated coalescent with recombination. Bioinformatics, 31(10), 1680–1682. https://doi.org/10.1093/bioinformatics/btu861


