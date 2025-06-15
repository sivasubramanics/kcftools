# `kcftools kcf2tsv`

The `kcftools kcf2tsv` command is used to convert a `.kcf` file into a `.tsv` file. The resulting output mimics the format used by IBSpy. This is typically used to export a single sample's data in a flat, tab-delimited format.

---

## Usage

    $ kcftools kcf2tsv -i=<kcfFile> -o=<outFile> [-s=<sampleName>]

---

## Options

| Option                         | Description                                                       | Required |
|--------------------------------|-------------------------------------------------------------------|----------|
| `-i`, `--input=<kcfFile>`      | Input KCF file name                                               | Yes      |
| `-o`, `--output=<outFile>`     | Output file name prefix                                           | Yes      |
| `-s`, `--sample=<sampleName>`  | Sample name to embed in the output (optional)                     | No       |


---

## Example

    $ kcftools kcf2tsv -i=data/sample.kcf -o=out/sample.tsv -s=S1

This command reads the `data/sample.kcf` file, converts it into a TSV file, and saves it to `out/sample.tsv`. The sample name `S1` is included in the output.

---
