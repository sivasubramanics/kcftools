<!--- badge: start --->
[![GitHub all releases](https://img.shields.io/github/downloads/sivasubramanics/kcftools/total?style=social&logo=github&label=Download)](https://github.com/sivasubramanics/kcftools/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/kcftools.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/kcftools)
[![Release](https://github.com/sivasubramanics/kcftools/actions/workflows/release.yml/badge.svg)](https://github.com/sivasubramanics/kcftools/actions/workflows/release.yml)
[![Version](https://img.shields.io/badge/version-0.0.1-green.svg)](https://github.com/sivasubramanics/kcftools/releases)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!--- badges: end --->

# KCFTOOLS

**KCFTOOLS** is a Java-based toolset for identifying genomic variations through counting kmer presence/absence between reference and query genomes. It utilizes precomputed k-mer count databases (from [KMC](https://github.com/refresh-bio/KMC)) to perform a wide array of genomic analyses including variant detection, IBS window identification, and genotype matrix generation.

---

## Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Limitation and Performance Notes](#-limitation-and-performance-notes)
- [Usage](#usage)
   - [getVariations](#getvariations)
   - [cohort](#cohort)
   - [findIBS](#findibs)
   - [splitKCF](#splitkcf)
   - [getAttributes](#getattributes)
   - [genotypeMatrix](#genotypematrix)
- [LICENSE](#license)
- [Contact](#contact)

---

## Introduction

KCFTOOLS is designed for high-throughput genomic analysis using efficient k-mer based methods. By leveraging fast k-mer counting from tools like KMC, KCFTOOLS can rapidly compare genome samples to a reference, identify variations, and produce downstream outputs useful for population genetics and comparative genomics studies.

---

## Features

- **Screen for Variations**: Detect sequence variations by comparing k-mers from reference and sample.
- **Cohort Creation**: Merge multiple `.kcf` sample files into a unified cohort.
- **IBS Window Identification**: Identify Identity-by-State (IBS) windows or variable regions across samples.
- **Chromosome-wise Splitting**: Partition KCF files by chromosome for parallel or targeted analysis.
- **Attribute Extraction**: Generate summaries and detailed statistics from `.kcf` files.
- **Genotype Matrix Generation**: Convert `.kcf` files into population-level genotype matrices.
- **Window Composition**: Compose larger genomic windows from finer-grained `.kcf` data.
- **Conversion Utilities**: Export `.kcf` files to TSV format (to replicate IBSpy output).

---

## Installation

You can install **kcftools** using either **Bioconda** or from source.

### 1. Using Bioconda (recommended) [currently under dev]

If you have Bioconda set up, simply run:

```bash
conda install -c bioconda kcftools
```

### 2. From Source

#### Requirements

- **Java 17+**
- **Maven** (for building)

#### Steps

1. Clone the repository:

   ```bash
   git clone https://github.com/sivasubramanics/kcftools.git
   cd kcftools
   ```

2. Build the project using Maven:

   ```bash
   mvn clean package
   ```

3. The JAR file will be located in the `target` directory:

   ```bash
   ls target/kcftools-<version>.jar
   ```

4. Run the tool:

   ```bash
   java -jar target/kcftools-<version>.jar <command> [options]
   ```

---

## ⚠️ Limitation and Performance Notes

1. **KMC DB Compatibility**:
    - `getVariations` plugin works only with KMC databases produced by `kmc` version 3.0.0 or higher. 
    - This version currently supports only KMC database files generated with a **signature length of 9** (i.e., using `-p 9`).  
   Files created with other signature lengths are not guaranteed to work and may lead to unexpected behavior.

2. **Memory Usage with `--memory` or `-m` Option**:  
   The `getVariations` plugin can be significantly faster when used with the `--memory` option, which loads the KMC database entirely into memory.  
   However, this may lead to Java heap space errors on large DBs. To prevent such issues:
    - Run with a custom heap size using the `-Xmx` JVM option  
      _Example: `kcftools -Xmx16G getVariations ...`_
    - Or, set the default heap size via the environment variable `KCFTOOLS_HEAP_SIZE`  
      _Example: `export KCFTOOLS_HEAP_SIZE=16G`_

---

## Usage

KCFTOOLS provides several subcommands. General usage:

```bash
kcftools <command> [options]
```

---
### `getVariations`
Detect and count variations by comparing reference k-mers with a query KMC database.

```bash
kcftools getVariations [options]
```

**Required Options:**

	-r, --reference=<refFasta>    : Reference FASTA file  
	-k, --kmc=<kmcDBprefix>       : KMC database prefix  
	-o, --output=<outFile>        : Output `.kcf` file  
	-s, --sample=<sampleName>     : Sample name  
	-f, --feature=<featureType>   : Feature granularity: `window`, `gene`, or `transcript`  

**Optional:**

	-t, --threads=<n>             : Number of threads (default: 2)  
	-w, --window=<size>           : Window size if `featureType=window`  
	-g, --gtf=<gtfFile>           : GTF annotation file (for gene/transcript features)  
	--wi, --wt, --wr              : Weights for inner distance, tail distance, and kmer ratio, respectively  
	-m, --memory                  : Load KMC database into memory (faster for small DBs)
---

### `cohort`
Combine multiple `.kcf` files into a single cohort for population-level analysis.

```bash
kcftools cohort [options]
```

**Options:**
   
      -i, --input=<file1>,<file2>,...  : Comma-separated list of KCF files
      -l, --list=<listFile>            : File containing newline-separated KCF paths
      -o, --output=<outFile>           : Output cohort `.kcf` file

---

### `findIBS`
Identify Identity-by-State (IBS) or variable regions in a sample.

```bash
kcftools findIBS [options]
```

**Options:**

      -i, --input=<kcfFile>      : Input KCF file
      -r, --reference=<refFasta> : Reference FASTA
      -o, --output=<outFile>     : Output `.kcf` file
      --bed                      : Also output BED file format
      --summary                  : Write summary TSV report
      --min=<minConsecutive>     : Minimum consecutive window count
      --score=<cutOff>           : Score threshold
      --var                      : Detect variable regions instead of IBS

---

### `splitKCF`
Split a KCF file by chromosome.

```bash
kcftools splitKCF [options]
```

**Options:**

      -k, --kcf=<kcfFile>     : Input KCF file
      -o, --output=<outDir>   : Output directory

---

### `getAttributes`
Extract and summarize metadata and attributes from a KCF file.

```bash
kcftools getAttributes [options]
```

**Options:**

      -i, --input=<kcfFile>   : Input `.kcf` file
      -o, --output=<prefix>   : Output prefix (produces `.tsv`, `.json`, etc.)

---

### `genotypeMatrix`
Generate a genotype matrix from a `.kcf` file, suitable for GWAS or population studies.

```bash
kcftools genotypeMatrix [options]
```

**Options:**
      
      -i, --input=<kcfFile>     : Input `.kcf` file
      -o, --output=<prefix>     : Output matrix prefix

---

## Notes

⚠️ **Warning:** Unit tests and comprehensive validation are still under development. Use results carefully and validate downstream.

---
## License

This project is licensed under the [GNU General Public License v3.0 only (GPL-3.0-only)](https://www.gnu.org/licenses/gpl-3.0.html).

See the [LICENSE](./LICENSE) file for details.

----

## Contact

For questions or contributions, please reach out to: c.s.sivasubramani@gmail.com

