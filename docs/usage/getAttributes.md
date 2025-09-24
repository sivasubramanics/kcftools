# `kcftools getAttributes`

The `getAttributes` command in `kcftools` is used to **extract attribute information** from a KCF (K-mer Comparison Format) file. This is useful for summarizing or analyzing the metadata or numerical properties computed during variation or IBS detection.

---

## Usage

    $ kcftools getAttributes -i=<kcfFile> -o=<outFile>

---

## Description

KCF files store per-window or per-feature variation or similarity scores along with associated metadata. `getAttributes` parses this data and outputs attribute summaries or tables in a more accessible format (e.g., TSV or similar) for downstream analysis.

The command takes an input `.kcf` file and produces one or more output files prefixed with the name you specify. These files may contain:

- Per-window attribute scores (e.g., variation score, IBS score, k-mer counts)
- Any other feature-level statistics computed during KCF generation

This can be used to perform statistical analyses, generate custom plots, or inspect data before further processing.

---

## Options

| Option                         | Description                                          | Required |
|--------------------------------|------------------------------------------------------|------|
| `-i`, `--input=<kcfFile>`      | Input KCF file from which to extract attributes      | Yes |
| `-o`, `--output=<outFile>`     | Output file **prefix** for saving extracted results  | Yes |
| `-a`, `--attributes=<attrList>` | Comma-separated list of specific attributes to extract (default: all) | No |

---

## Output

The output will typically include one or more `.tsv` or `.txt` files prefixed by the provided output name. The actual file(s) depends on the structure of the input KCF and may contain:

- Region identifiers (chromosome, start, end)
- Per-sample or per-region attribute scores
- Metadata such as sample name, feature type, score thresholds

For example:

    $ attributes_output.tsv

---

## Example

    $ kcftools getAttributes -i sample.kcf -o sample_attributes

This will extract all available attribute data from `sample.kcf` and write it to files prefixed with `sample_attributes`.

---

!!! note
    - This command is **read-only**—it does not modify the original KCF.
    - Useful when preparing data for visualization (e.g., in R or Python plotting tools).
    - You can use this to summarize scores before applying filters with tools like `findIBS` or `kcfToMatrix`.

---

## Help

To view help for this command:

    $ kcftools getAttributes --help
