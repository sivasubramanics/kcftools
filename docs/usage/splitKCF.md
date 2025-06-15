# `kcftools splitKCF`

The `splitKCF` command in `kcftools` is used to **split a single KCF file into per-chromosome files**. This is especially useful when working with large genome-wide variation data and you need to analyze or visualize data one chromosome at a time.

---

## Usage

    $ kcftools splitKCF -k=<kcfFile> -o=<outDir>

---

## Description

KCF (K-mer Comparison Format) files may contain genome-wide variation or similarity data. When working with such large-scale data, it is often necessary to analyze data chromosome by chromosome, or to parallelize processing by splitting the input per chromosome.

This command reads a full KCF file and generates a **separate KCF file for each chromosome** found in the original input. The split files are saved in the specified output directory, and named appropriately (usually including the chromosome name or number).

---

## Options

| Option                        | Description                                 | Required |
|-------------------------------|---------------------------------------------|----------|
| `-k`, `--kcf=<kcfFile>`       | Input KCF file to be split                  | Yes      |
| `-o`, `--output=<outDir>`     | Output directory where files will be saved  | Yes      |

---

## Output

For each chromosome in the input KCF file, a new `.kcf` file will be created in the specified output directory. For example, if your input file contains data for chromosomes `chr1`, `chr2`, and `chrX`, the output directory will contain:

    output_directory/
    ├── chr1.kcf
    ├── chr2.kcf
    └── chrN.kcf

These files can then be used for chromosome-specific analysis, visualization, or downstream pipelines such as matrix conversion (`kcfToMatrix`), plotting, or region filtering.

---

## Example

    $ kcftools splitKCF -k merged_variations.kcf -o split_by_chrom/

This command will read `merged_variations.kcf` and generate one file per chromosome under the `split_by_chrom/` directory.

---

!!! note
    - Ensure the input KCF file is correctly formatted and includes chromosome information per entry.
    - The output directory will be created if it does not already exist.
    - This is a utility command—useful before plotting, matrix generation, or when managing large data.

---

## Help

To view help for this command:

    $ kcftools splitKCF --help
