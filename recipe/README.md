# Locally building via `conda-build`

Use this recipe for locally creating your own `conda` package. 

## Operating Systems:

- **Linux**: Tested, verified.
- **MacOSX**: Work in progress, but follow [this](https://docs.conda.io/projects/conda-build/en/latest/resources/compiler-tools.html#macos-sdk) guide for installing the SDK which is *not* packaged with `conda-build`. 
- **Windows**: Currently no support, nor is it planned.

## Usage

To build your package (after cloning the repository and installing `conda-build`)

```sh
conda build . -c terhorst
```

To install:

```sh
conda install -c ${CONDA_PREFIX}/conda-bld smcsmc
```

This will install both the python package and `smcsmc`/`scrm` binaries. 
