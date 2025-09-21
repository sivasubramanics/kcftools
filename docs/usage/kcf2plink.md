# kcftools kcf2plink

The `kcftools kcf2plink` command is used to convert a `.kcf` file into PLINK format files (`.ped` and `.map`). This allows users to leverage the KCF data for further genetic analyses using PLINK.

---

## Usage

    $ kcftools kcf2plink [-a=<scoreA>] [-b=<scoreB>] [--chrs=<chrsFile>]
                          -i=<inFile> [--maf=<minMAF>]
                          [--max-missing=<maxMissing>] -o=<outPrefix>
                          [--score_n=<scoreN>] 
---
## Description

The `kcf2plink` command processes a KCF file to generate PLINK-compatible files. It interprets the identity scores in the KCF file to classify genomic regions into reference, alternate, or missing categories based on user-defined thresholds. The output files can then be used for various genetic analyses, including association studies and population genetics.

## Options
| Option                                | Description                                                                 | Default    | Required |
|---------------------------------------|-----------------------------------------------------------------------------|------------|----------|
| `-i`, `--input=<inFile>`              | Input `.kcf`
    file to convert                                                | _N/A_      | Yes      |
| `-o`, `--output=<outPrefix>`          | Prefix for output files (PED and MAP)                                     | _N/A_      | Yes      |
| `-a`, `--score_a=<scoreA>`            | Lower score cut-off
    to call the **reference allele**                        | *none*     | No       |
| `-b`, `--score_b=<scoreB>`            | Lower score cut-off
    to call the **alternate allele**                        | *none*     | No       |
| `-n`, `--score_n=<scoreN>`            | Upper score cut-off
    to call a region as **missing** (between `score_b` and `score_n`) | *none*     | No       |
| `--maf=<minMAF>`                      | Minimum minor allele frequency to retain a region                           | *none*     | No       |
| `--max-missing=<maxMissing>`          | Maximum proportion of missing data allowed per region                       | *none*     | No       |
| `--chrs=<chrsFile>`                   | File listing chromosomes to include (one per line)                          | *all*      | No       |

---
## Output

This feature generates two output files:
- A PED file (`<outPrefix>.ped`) containing genotype information for each sample.
- A MAP file (`<outPrefix>.map`) containing marker information.
- These files are formatted according to PLINK specifications and can be directly used in PLINK for further analyses.

**Example output files:**

    output.ped
    output.map

---
## Example

**Convert a KCF to PLINK format with thresholds:**

    $ kcftools kcf2plink -i sample.kcf -o plink_output -a 95 -b 70 --maf 0.05 --max-missing 0.1
This command:
- Treats windows with score ≥95 as **reference**
- Treats windows with score ≥70 (but <95) as **alternate**
- Filters out regions with MAF <5% or >10% missing data
- Generates `plink_output.ped` and `plink_output.map` files for use in PLINK

---

!!! note
    - **Score thresholds (`--score_a`, `--score_b`)** are crucial to define genotype states.
    - **MAF** and **missing data** thresholds help clean noisy or non-informative regions.
    - Use the `--chrs` option to restrict to a subset of chromosomes if desired.
    - If allele scores are ambiguous or missing for a sample in a region, that entry will be marked `-1`.

---

## Help

To get help on this command, run:

    $ kcftools kcf2plink --help