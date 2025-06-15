# `kcftools findIBS`

The `findIBS` command in `kcftools` is used to identify **Identical-By-State (IBS)** or **Variable Regions** in a KCF (K-mer Comparison Format) file. The type of region detected depends on whether the `--var` flag is used.

---

## Usage

    $ kcftools findIBS [--bed] [--summary] [--var] -i=<inFile> 
                 [--min=<minConsecutive>] -o=<outFile> 
                 [--score=<scoreCutOff>]

---

## Options

| Option                     | Description                                                                 | Default       |
|---------------------------|-----------------------------------------------------------------------------|---------------|
| `-i`, `--input=<inFile>`  | Input KCF file name (required)                                              | _N/A_         |
| `-o`, `--output=<outFile>`| Output KCF file name (required)                                             | _N/A_         |
| `--min=<minConsecutive>`  | Minimum number of consecutive windows required for detection                | `4`           |
| `--score=<scoreCutOff>`   | Score cut-off percentage for defining IBS/variable windows                  | `95.00`       |
| `--bed`                   | Write output in BED file format                                             | `false`       |
| `--summary`               | Write summary TSV file                                                      | `false`       |
| `--var`                   | Detect **Variable Regions** instead of IBS regions                          | `false`       |

---

## Description

By default, `kcftools findIBS` detects **IBS windows**, which are stretches of the genome where samples are highly similar (based on the score threshold). When the `--var` flag is used, it instead detects **Variable Regions**, where the similarity score falls *below* the threshold—indicating higher variation.

### IBS vs Variable Regions

| Mode         | Meaning                                                                                   |
|--------------|--------------------------------------------------------------------------------------------|
| **IBS**      | Detects windows with similarity **above** the threshold, indicating shared genetic regions |
| **Variable** | Detects windows with similarity **below** the threshold, indicating higher divergence      |

---

## Example

    $ kcftools findIBS --input=input.kcf --output=ibs_regions.kcf --score=98.0 --min=5

This command finds IBS regions where the similarity score is ≥98.0 and stretches are at least 5 consecutive windows long.

    $ kcftools findIBS --input=input.kcf --output=variable_regions.kcf --var --score=90.0

This command finds variable regions (regions with similarity <90.0).

---

## Output

Depending on the options used, `kcftools findIBS` may produce:

- A filtered `.kcf` file (`--output`)
- A BED format file (`--bed`)
- A summary `.tsv` file (`--summary`)
