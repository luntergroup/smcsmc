# Locally building via `conda-build`

Use this recipe for locally creating your own `conda` package. It only works on Linux right now.

## Usage

To build your package, you need quite a few dependencies specific to building conda packages. I've placed my (working) build environment in `build.yml`, which you can import to your own system.  

Once this has been succesful, on Linux you can build the package with: 

```sh
conda build .
```

To install:

```sh
conda install -c ${CONDA_PREFIX}/conda-bld smcsmc
```

This will install both the python package and `smcsmc`/`scrm` binaries. 
