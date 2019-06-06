# SMCSMC
## Demographic Inference using Particle Filters

SMCSMC (Sequential Monte Carlo for the Sequential Markovian Coalescent) or SMC2 is a program for inferring population history from multiple genome sequences. 

## Installation

This repository contains two effectively independant components, and both must be installed to properly use SMC2. We have automated this process in a `conda` package, and we highly recommend installing it this way.

```sh
conda install -c terhorst, conda-forge smcsmc
```

### Installation from Source

Alternatively, a combination of `cmake` and `pip` can be used to install the python and core components.

#### To install the backend

```sh
mkdir build; cd build
cmake ..
cmake
```

#### To install the frontend

```sh
pip install .
```
