# Changelog

This file contains notable updates to `smcsmc` partitioned by version.

## [1.1.0](https://github.com/luntergroup/smcsmc/compare/1.0.2...v1.1.0) (2022-10-18)


### Features

* basic docker image, no data ([665ddbf](https://github.com/luntergroup/smcsmc/commit/665ddbf275bf36420b5fe8c604093a1cd1ecd5ed))
* build smcsmc on circle  ([#75](https://github.com/luntergroup/smcsmc/issues/75)) ([3efa55e](https://github.com/luntergroup/smcsmc/commit/3efa55ece17a2725add518bbc3bf1e0c46f1a832))
* circle badge and release please action ([66f118d](https://github.com/luntergroup/smcsmc/commit/66f118d9e1a3cd656545b345828411032849f87c))
* test cli and add circle for python package ([90c5238](https://github.com/luntergroup/smcsmc/commit/90c5238eb8796de3f173afdf521dbad91a62922f))


### Bug Fixes

* move cli logic to an entry point and check deps ([fd94f9e](https://github.com/luntergroup/smcsmc/commit/fd94f9e1063c6df63e4b5a87d9f122c74baace0e))


### Documentation

* docker documentation ([14b7325](https://github.com/luntergroup/smcsmc/commit/14b7325968f6e0dbaafe064a28d5e0217369681c))

## [Unreleased]

- Fixed a bug where the `scrm` approximation window was not being properly set.

## [1.0.2] - 16/10/19
### Added
- Scripts to convert from the ARGs output by `smcsmc` to a partial `tskit` tree-sequence. Notably, this is used to infer better migration tracts.
- Build environment for local `conda` building.

### Changed
- Default parameters of simulations, in order to correctly bound the proportions to the correct epochs.
- Format of the simulations to be more expressive.

### Fixed
- `rpath` of `smcsmc` when installed now automatically prefers `conda` libraries.
