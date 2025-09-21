# Changelog

All notable changes to this project will be documented in this file.

---

## [v0.3.0] - 2025-09-21

### Added
- `kcf2gt` plugin to convert KCF files to genotype table format for tools like Tassel and GAPIT [`#2c2c4a3`](https://github.com/sivasubramanics/kcftools/commit/2c2c4a3b84ae69d48b8370cf3416ae2671917d91)
- `scoreRecalc` plugin to recalculate scores in KCF files based on new parameters [`#8f60682`](https://github.com/sivasubramanics/kcftools/commit/8f606826530cf3210e8ad502594f54e55f2f6bc5)
- Addes sliding window option to `getVariations` [`#4f66d77`](https://github.com/sivasubramanics/kcftools/commit/4f66d77350b6b8e4f83fe52f0dadf08052f1d9d9)

### Fixed
- Fixed hard coded score_a and score_b values in `kcf2gt` nad `kcfToPed` plugins [`#1bdddfd`](https://github.com/sivasubramanics/kcftools/commit/1bdddfdce8feecb4ea6a8378fe5a8650e37e0f84)
- Fixed min consecutive windows condition in `findIBS` [`#41f5ecd`](https://github.com/sivasubramanics/kcftools/commit/41f5ecd3371e64db848388ea072cddb34da8aeb5)
- Fixed weights handling in `getVariations` [`#76d8d79`](https://github.com/sivasubramanics/kcftools/commit/76d8d796c01cc088dad5fb0d17c00a7dd83d0419)

### Changed
- Removed plugin `kcfToMatrix` in favor of `kcf2gt` [`#1bdddfd`](https://github.com/sivasubramanics/kcftools/commit/1bdddfdce8feecb4ea6a8378fe5a8650e37e0f84)
- Disabled `CompareIBS` plugin [`#e3a747a`](https://github.com/sivasubramanics/kcftools/commit/e3a747a48058dc0c0a45a06c4f8271ef4bf16b9c)


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

