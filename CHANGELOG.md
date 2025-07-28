# Changelog

All notable changes to this project will be documented in this file.

---

## [v0.2.0] - 2025-07-28
### Added
- `kcfToPed` experimental plugin to convert KCF files to PLINK PED format [`#37f1729`](https://github.com/sivasubramanics/kcftools/commit/37f1729f973f0814c8724d9f2dfd0693cbc8ebf5)
- Added readthedocs documentation [`#5bf98a1`](https://github.com/sivasubramanics/kcftools/commit/5bf98a1787e5f2ba9232872a20e1d78cb5f52a14)

### Fixed
- Fixed tail distance and total_kmer calculation in `increaseWindows` [`#4a76768`](https://github.com/sivasubramanics/kcftools/commit/4a767681994ef26ea67204d1794c564be5864c80)

### Changed
- Updated GTF parsing in `getVariations` to handle AGAT formatted GTF files. [`#a68beb5`](https://github.com/sivasubramanics/kcftools/commit/a68beb59996880ceec31660a5db6c7c33cf35c32)

## [v0.1.0] - 2025-05-31
### Added
- `--maf` and `--maxmissing` parameters to `kcfToMatrix` [`#e3e3b90`](https://github.com/sivasubramanics/kcftools/commit/e3e3b90)
- `min-k-count` count parameter to `getVariations` [`#7aad1c1`](https://github.com/sivasubramanics/kcftools/commit/7aad1c1)
- Workflow and Methodology sections in README.md [`#fcf6c34`](https://github.com/sivasubramanics/kcftools/commit/fcf6c34), [`#8c7f83d`](https://github.com/sivasubramanics/kcftools/commit/8c7f83d)

### Fixed
- Fixed thread blocking issue in `getVariations` [`#7aad1c1`](https://github.com/sivasubramanics/kcftools/commit/7aad1c188f927b392d78df806524b351bd56d888) 

### Changed
- Updated documentation in README.md [`#c121ab2`](https://gitbub.com/sivasubramanics/kcftools/commit/c121ab2)

---

## [v0.0.1] - 2025-05-22
### Added
- First minor release with core functionality. [`#e465f4a`](https://github.com/sivasubramanics/kcftools/commit/e465f4a)
- Initial public release of `kcftools` in Bioconda.

