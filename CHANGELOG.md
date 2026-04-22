# Changelog

All notable changes to this project will be documented in this file.

---

## [2026-04-22]

### Added
- plugin to get Kmer density from multi fasta file [`#45aa140`](https://github.com/sivasubramanics/kcftools/commit/45aa140651bf4ee66850c1d1fac66dc596388ac2)
- added a attribute dist to return the tail + in dist [`#69a5c88`](https://github.com/sivasubramanics/kcftools/commit/69a5c8851c18810c4bb920e2bfe70e67a23ab0a1)
- script to convert kcf gt data into plink.ped format [`#2348160`](https://github.com/sivasubramanics/kcftools/commit/2348160b42e5a37b4566ba3e4d02e03be1abd576)
- example dataset [`#66e151e`](https://github.com/sivasubramanics/kcftools/commit/66e151e968fe3a0cd999126b9dcdc0c7e8a12c46)
- utility script to parse gtf [`#2f799e7`](https://github.com/sivasubramanics/kcftools/commit/2f799e7dd7efa66e496986d7ac70708034175f6d)
- meta information for the run findIBS in summary.tsv [`#c012d8f`](https://github.com/sivasubramanics/kcftools/commit/c012d8f5f25636b780a2fed5bf2e7f2028b7e44c)
- extractKCF feature [`#c44593f`](https://github.com/sivasubramanics/kcftools/commit/c44593fa9f15021f87c4b1d0fa91f775cfadc236)
- kmerCount or kmerDepth to the KCF file [`#1c32378`](https://github.com/sivasubramanics/kcftools/commit/1c32378d9461e6e65c72a58bb4a3c9674ad562f6)
- step based increaseWindow feature [`#6fdecbc`](https://github.com/sivasubramanics/kcftools/commit/6fdecbc73e325cbc3df6b0d09682c264169f7dba)
- multithreading support for cohorting [`#26295ad`](https://github.com/sivasubramanics/kcftools/commit/26295adfaabbd088ead9b6e108e6f8b30360a80d)
- file validation checkpoints [`#3c23a93`](https://github.com/sivasubramanics/kcftools/commit/3c23a9310494a7fd2d0ddc6ba9f2a60ce6147610)

### Fixed
- sentence about the contents [`#9449cbc`](https://github.com/sivasubramanics/kcftools/commit/9449cbc5db298762897c5f28ae20b0c51b6ea520)
- if older KCF version is used KD to be fixed [`#8feb203`](https://github.com/sivasubramanics/kcftools/commit/8feb20320cb885a58e91aa4da88e351b863e403d)
- code refactoring [`#5a18e24`](https://github.com/sivasubramanics/kcftools/commit/5a18e24f4f76919bf041d9434454e641aff6ce01)

### Changed
- KCFTOOLS to KCFtools in docs [`#cec5159`](https://github.com/sivasubramanics/kcftools/commit/cec51590296c35e2471a672a45bf58382ca12cb2)
- KCFTOOLS to KCFtools [`#90ed6ee`](https://github.com/sivasubramanics/kcftools/commit/90ed6ee8533e7fd0aa514f2525ca90e6d3c76352)
- added citation and more [`#07bcfb2`](https://github.com/sivasubramanics/kcftools/commit/07bcfb28658c4b6591ef8427f46fbd8cd8225f3b)
- now we can filter the IBS summary based on scores to plot [`#9d3ac44`](https://github.com/sivasubramanics/kcftools/commit/9d3ac442f8f3838cc86a5e03e28a2679cceee405)
- Added plugin to the main menu for Kmer density from multi fasta file [`#2246ed3`](https://github.com/sivasubramanics/kcftools/commit/2246ed39c84e96063e2b2003204952691586cab3)
- now we can choose the chrom based on the min length [`#b37c9ab`](https://github.com/sivasubramanics/kcftools/commit/b37c9abeef8fdd497109aefd4e144f5af58865e5)
- writing findIBS output in kcf is optional now [`#43d1c30`](https://github.com/sivasubramanics/kcftools/commit/43d1c30e0908d98a362f6ecb40b3b626f4380e59)
- removed unnecessary libraries [`#90724b6`](https://github.com/sivasubramanics/kcftools/commit/90724b6e20f68964f7aab8b3247678cf6c3da43b)
- added tree based ordering in the plotIBS.R [`#ac64b87`](https://github.com/sivasubramanics/kcftools/commit/ac64b87724bdf1bc8ac5af07297f1a4cb4c60518)
- move back to non multi threaded cohort [`#da627bc`](https://github.com/sivasubramanics/kcftools/commit/da627bc9509c5caa29b64341fcdd08ce06df9c5d)


## [v0.5.0] - 2025-12-20
### Added
 - Added example datasets for the pipeline in the `examples` directory [`#66e151e`](https://github.com/sivasubramanics/kcftools/commit/66e151e968fe3a0cd999126b9dcdc0c7e8a12c46)
 - Added utility script `GTF_utils.py` for GTF file manipulation and annotation of IBS summary files.[`#2f799e7`](https://github.com/sivasubramanics/kcftools/commit/2f799e7dd7efa66e496986d7ac70708034175f6d)
 - Added ``


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


## [v0.0.1] - 2025-05-22
### Added
- First minor release with core functionality. [`#e465f4a`](https://github.com/sivasubramanics/kcftools/commit/e465f4a)
- Initial public release of `kcftools` in Bioconda.
