# `kcftools kcf2gt`

The `kcf2gt` command in `kcftools` converts a **KCF (K-mer Comparison Format)** file into a **genotype table**, enabling downstream population genetics, phylogenetics, and statistical analysis. The resulting matrix encodes the allele state for each sample across windows or features, which can be further filtered or visualized.

---

## Usage

    $ kcftools kcf2gt [-r] [--score_a=<scoreA>] [--score_b=<scoreB>][--score_n=<scoreN>]
                     [--chrs=<chrsFile>] -i=<inFile>
                     [--maf=<minMAF>] [--max-missing=<maxMissing>]
                     -o=<outFile>

---

## Description

This tool transforms the region-level variation or IBS data from a `.kcf` file into a **genotype table**, where rows correspond to regions (windows, genes, etc.) and columns represent samples.

Each table entry is inferred based on scoring thresholds for the reference and alternate alleles, which are defined using `--score_a` and `--score_b`. You can also apply **variant filtering** using thresholds for **minor allele frequency (MAF)** and **missing data**.

The result is a matrix file (TSV) that can be used for:

- Principal component analysis (PCA)
- Phylogenetic trees
- Clustering and population structure analysis
- Association studies (e.g., GWAS)

---

## Options

| Option                                | Description                                                                 | Default    | Required |
|---------------------------------------|-----------------------------------------------------------------------------|------------|----------|
| `-i`, `--input=<inFile>`              | Input `.kcf` file to convert                                                | _N/A_      | Yes      |
| `-o`, `--output=<outPrefix>`          | Prefix for output files (matrix and metadata)                              | _N/A_      | Yes      |
| `-a`, `--score_a=<scoreA>`            | Lower score cut-off to call the **reference allele**                        | *none*     | No       |
| `-b`, `--score_b=<scoreB>`            | Lower score cut-off to call the **alternate allele**                        | *none*     | No       |
| `-n`, `--score_n=<scoreN>`            | Upper score cut-off to call a region as **missing** (between `score_b` and `score_n`) | *none*     | No       |
| `--maf=<minMAF>`                      | Minimum minor allele frequency to retain a region                           | *none*     | No       |
| `--max-missing=<maxMissing>`          | Maximum proportion of missing data allowed per region                       | *none*     | No       |
| `--chrs=<chrsFile>`                   | File listing chromosomes to include (one per line)                          | *all*      | No       |

---

## Output

Depending on options used, the command produces:

- A genotype table file (e.g., `.tsv`) containing region identifiers in row and sample data on columns
- Contig map file with numeric ID and matching chromosome/contig names

**Example output files:**

    output.tsv
    output.tsv.contigsMap.tsv

In the genotype table:

- Values typically represent allele status: `0` (reference), `2` (alternate), `-1` (missing), `1` (heterozygous, if applicable)
- Columns = sample names
- Rows = regions (windows, genes, etc.)

---

## Example

**Convert a KCF to genotype matrix with thresholds:**

    $ kcftools kcf2gt -i sample.kcf -o matrix_output -a 95 -b 70 --maf 0.05 --max-missing 0.1

This command:

- Treats windows with score ≥95 as **reference**
- Treats windows with score ≥70 (but <95) as **alternate**
- Filters out regions with MAF <5% or >10% missing data

---

!!! note
    - **Score thresholds (`--score_a`, `--score_b`)** are crucial to define genotype states.
    - **MAF** and **missing data** thresholds help clean noisy or non-informative regions.
    - Use the `--chrs` option to restrict to a subset of chromosomes if desired.
    - If allele scores are ambiguous or missing for a sample in a region, that entry will be marked `-1`.

---

## Help

To view the command-line help:

    $ kcftools kcf2gt --help
