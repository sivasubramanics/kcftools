# `kcftools increaseWindow`

The `increaseWindow` command in `kcftools` is used to **aggregate adjacent windows** in a `.kcf` file into **larger non-overlapping windows** of a specified size. This operation can help smooth the signal, reduce noise, or adapt the resolution of KCF-based variation or similarity analyses to broader genomic scales.

---

## Usage

    $ kcftools increaseWindow -i=<inFile> -o=<outFile> -w=<windowSize>

---

## Description

KCF files often contain fine-scale data, typically with window sizes like 100 bp or 1 kb. However, certain analyses—such as population-level statistics, whole-chromosome summaries, or broad variation scans—require coarser resolution.

This command **merges consecutive windows** to form **larger, aggregated windows** of the user-specified size. For example, if the original KCF file uses 1,000 bp windows, and `--window 5000` is specified, every 5 adjacent windows (on the same chromosome) will be grouped into one.

During merging, attributes like variation scores, IBS values, and *k*-mer counts are **aggregated**, typically by recomputing, depending on the attribute type.

---

## Options

| Option                          | Description                                                | Required |
|---------------------------------|------------------------------------------------------------|-----|
| `-i`, `--input=<inFile>`        | Input `.kcf` file to process                               | Yes |
| `-o`, `--output=<outFile>`      | Output KCF file with merged windows                        | Yes |
| `-w`, `--window=<windowSize>`   | Desired new window size (in base pairs)                    | Yes |

---

## Output

The output is a `.kcf` file containing fewer but larger windows. Each new window is a merged region that summarizes attributes from its constituent smaller windows.

---

## Example

    $ kcftools increaseWindow -i genome_1kb.kcf -o genome_10kb.kcf -w 10000

This merges the original 1,000 bp windows into 10,000 bp windows, reducing granularity and improving interpretability at broader scales.

---

!!! note
    - Windows are merged **only if they are contiguous and on the same chromosome**.
    - The final window may be smaller if the chromosome length is not a perfect multiple of the new window size.
    - This operation is ideal **after variation detection** (e.g., from `getVariations`) and **before matrix conversion** (`kcfToMatrix`) when coarser analysis is desired.

---

## Help

To view command-line help:

    $ kcftools increaseWindow --help
