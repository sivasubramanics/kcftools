# Welcome to KCFTOOLS Documentation!

**KCFTOOLS** is a Java-based toolset for identifying genomic variations through counting kmer presence/absence between reference and query genomes. It utilizes precomputed *k*-mer count databases (from [KMC](https://github.com/refresh-bio/KMC)) to perform a wide array of genomic analyses including variant detection, IBS window identification, and genotype matrix generation.

---
## License

This project is licensed under the [GNU General Public License v3.0 only (GPL-3.0-only)](https://www.gnu.org/licenses/gpl-3.0.html).

---

## Contents
### General Information
- [Quick Start](general/quickstart.md)
- [Installation](general/installation.md)
- [Workflow](general/workflow.md)
- [Features](general/features.md)
- [Methodology](general/methodology.md)
### Usage
- [Command Line Interface](usage/cli.md)
- [getVariations](usage/getVariations.md)
- [cohort](usage/cohort.md)
- [findIBS](usage/findIBS.md)
- [splitKCF](usage/splitKCF.md)
- [getAttributes](usage/getAttributes.md)
- [kcf2gt](usage/kcf2gt.md)
- [kcf2tsv](usage/kcf2tsv.md)
- [increaseWindow](usage/increaseWindow.md)
- [scoreRecalc](usage/scoreRecalc.md)
- [kcf2plink](usage/kcf2plink.md)
### File Formats
- [KCF File](formats/kcf.md)
- [KMC Database](formats/kmc.md)
- [Genotype Table File](formats/gttable.md)
- [Attributes File](formats/attributes.md)
- [IBS Summary File](formats/ibs.md)
### [Additional Information](misc.md)