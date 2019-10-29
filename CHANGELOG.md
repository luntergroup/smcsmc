# Changelog

This file contains notable updates to `smcsmc` partitioned by version.

## [Unreleased]

## [1.0.2] - 16/10/19
### Added
- Scripts to convert from the ARGs output by `smcsmc` to a partial `tskit` tree-sequence. Notably, this is used to infer better migration tracts.
- Build environment for local `conda` building.

### Changed
- Default parameters of simulations, in order to correctly bound the proportions to the correct epochs.
- Format of the simulations to be more expressive.

### Fixed
- `rpath` of `smcsmc` when installed now automatically prefers `conda` libraries. 
