# `kcftools getVariations`

The `getVariations` command screens for reference k-mers that are **not** present in the KMC database, identifying candidate variations.

---

## Usage

    $ kcftools getVariations [OPTIONS]

---

## Required Options

| Option                           | Description                                                                 |
|----------------------------------|-----------------------------------------------------------------------------|
| `-r`, `--reference=<refFasta>`   | Reference FASTA file                                                        |
| `-k`, `--kmc=<kmcDBprefix>`      | KMC database prefix (omit `.kmc_pre` and `.kmc_suf`)                        |
| `-o`, `--output=<outFile>`       | Output file name (in KCF format)                                            |
| `-s`, `--sample=<sampleName>`    | Sample name to associate with the output                                    |
| `-f`, `--feature=<featureType>`  | Feature type for variation detection: `window`, `gene`, or `transcript`    |

---

## Optional Parameters

| Option                                 | Description                                                                 | Default  |
|----------------------------------------|-----------------------------------------------------------------------------|----------|
| `-t`, `--threads=<nThreads>`           | Number of threads to use                                                    | `2`      |
| `-m`, `--memory`                       | Load entire KMC DB into memory for faster access                            | `false`  |
| `--wi=<innerDistanceWeight>`           | Weight for inner *k*-mer distance in scoring                                | `0.3`    |
| `--wt=<tailDistanceWeight>`            | Weight for tail *k*-mer distance in scoring                                 | `0.3`    |
| `--wr=<kmerRatioWeight>`               | Weight for *k*-mer ratio in scoring                                         | `0.4`    |
| `-w`, `--window=<windowSize>`          | Window size in base pairs (used with `--feature=window`)                    | _N/A_    |
| `-g`, `--gtf=<gtfFile>`                | GTF file with annotations (required for `gene` or `transcript` features)   | _N/A_    |
| `-c`, `--min-k-count=<minKmerCount>`   | Minimum *k*-mer count threshold to consider valid                           | `1`      |
| `-p`, `--step=<stepSize>`                | Step size in base pairs for sliding windows (used with `--feature=window`)  | `windowSize` |

---

## Examples

**Detect variation in 1000 bp windows:**
  
    $ kcftools getVariations -r ref.fa -k sample_kmc -o sample.kcf -s sample_name -f window -w 1000

**Detect variation in genes using GTF:**
    
    $ kcftools getVariations -r ref.fa -k sample_kmc -o sample.kcf -s sample_name -f gene -g annotations.gtf

**Custom *k*-mer weights and multithreading:**

    $ kcftools getVariations -r ref.fa -k sample_kmc -o sample.kcf -s sample_name -f window -w 1000 --wr 0.5 --wi 0.2 --wt 0.3 -t 8 -m

---

!!! note
- Ensure your KMC database was generated with compatible parameters (e.g., *k*-mer size) relative to the reference genome.
- Use appropriate `--feature` settings depending on your annotation and biological context.
- The output `.kcf` file can be used for downstream analysis with other `kcftools` commands like `findIBS` or `kcfToMatrix`.

---

## Help

To display command-line help:

    $ kcftools getVariations --help
