# Locally building via `conda-build`

Use this recipe for locally creating your own `conda` package. It only works on Linux right now.

## Usage

You can build the package with 

```sh
conda build .
```

To install from the local build files:

```sh
conda install -c conda-forge -c ${CONDA_PREFIX}/conda-bld smcsmc
```

This will install both the python package and `smcsmc`/`scrm` binaries. 

## Building for other Python versions

Build variants are included in the `conda_build_config.yaml`. The only real constraint on the python versions that we are able to build for is the `boost` package. Note that the version on `conda-forge` is more current than the default channels, so this should be set preferentially on the command line.
