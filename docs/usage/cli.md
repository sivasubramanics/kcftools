# `kcftools` Command Line Interface

`kcftools` provides a suite of tools to handle k-mer counting based variation analysis.

## Usage

    $ kcftools [-hV] [COMMAND]

### Options

- `-h`, `--help` : Show this help message and exit.
- `-V`, `--version` : Print version information and exit.

---

## Available Commands

### `getVariations`

Screen for reference kmers that are not present in the KMC database, and detect variation.

    $ kcftools getVariations [OPTIONS]

### `cohort`

Create a cohort of sample KCF files.

    $ kcftools cohort [OPTIONS]

### `findIBS`

Find Identity-By-State (IBS) windows in a KCF file.

    $ kcftools findIBS [OPTIONS]

### `splitKCF`

Split a KCF file by chromosome.

    $ kcftools splitKCF [OPTIONS]

### `getAttributes`

Extract attributes from KCF files.

    $ kcftools getAttributes [OPTIONS]

### `kcfToMatrix`

Convert KCF windows into a genotype matrix.

    $ kcftools kcfToMatrix [OPTIONS]

### `kcf2tsv`

Convert a KCF file into a TSV file (similar to IBSpy output).

    $ kcftools kcf2tsv [OPTIONS]

### `increaseWindow`

Increase the window size of a KCF file by merging windows.

    $ kcftools increaseWindow [OPTIONS]

---

_For detailed help on a specific command, run:_

    $ kcftools <COMMAND> --help
