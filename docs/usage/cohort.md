# `kcftools cohort`

The `cohort` command allows you to **combine multiple KCF files** (each representing an individual sample) into a **single cohort file**. This is essential for performing joint analyses across multiple samples, such as identifying shared or unique variation patterns.

## Usage

    $ kcftools cohort [OPTIONS]

## Required Option

- `-o`, `--output=<outFile>`  
  Output file name for the resulting cohort KCF file.

## Input Options (one of the following is required)

You must specify either `--input` or `--list` to provide sample KCF files.

- `-i`, `--input=<inFiles>[,<inFiles>...]`  
  A comma-separated list of individual sample `.kcf` files.  
  Example: `sample1.kcf,sample2.kcf,sample3.kcf`

- `-l`, `--list=<listFile>`  
  A plain text file with one `.kcf` file path per line.

---

## Example Usages

### Combine samples using `--input`

    $ kcftools cohort \
    --input sample1.kcf,sample2.kcf,sample3.kcf \
    --output cohort.kcf

This command will merge the three sample `.kcf` files into one file named `cohort.kcf`.

---

### Combine samples using a list file

Suppose you have a file `sample_list.txt` with the following contents:

```
sample1.kcf
sample2.kcf
sample3.kcf
```

Run the following:

    $ kcftools cohort \
    --list sample_list.txt \
    --output cohort.kcf

This is especially useful when working with large cohorts of samples.

---

!!! note
    - All input `.kcf` files should be generated using the **same reference**, **k-mer size**, and **feature type** to ensure compatibility.
    - The output cohort file (`cohort.kcf`) can be used in downstream analyses such as IBS calculation with `kcftools findIBS` or genotype matrix conversion with `kcftools kcfToMatrix`.
    - You cannot use both `--input` and `--list` at the same time — choose one method for specifying input files.
