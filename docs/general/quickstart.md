# Quick Start

This is a quick start guide to use **KCFTOOLS** for genomic variation analysis. It covers installation, database preparation, and basic usage. We recommend reading the [Workflow](/general/workflow/) and [Methodology](/general/methodology/) sections for a more detailed understanding of the tool's functionalities.

---

## 1. Install KCFTOOLS

Assuming the prerequisites are met, you can install KCFTOOLS using the following command:

    $ conda install -c bioconda kcftools

## 2. Prepare KMC Database

Ensure you have a KMC database ready. If you don't have one, you can create it using KMC3:

    $ kmc -k31 -ci1 -cs10000 -fm -t4 @reads.txt sample1_kmc tmp/

- @reads.txt contains a list of input FASTQ files.
- Output will be sample1_kmc.kmc_pre and sample1_kmc.kmc_suf.

## 3. Create variations KCF file
To create a KCF file containing variations, use the `getVariations` command:

    $ kcftools getVariations \
        --kmc sample.kmc_db \
        --reference reference.fasta \
        --sample sample \
        --output sample.kcf \
        --window 10000 \
        --feature window \
        --memory \
        --threads 2

This command will generate a KCF file named `sample.kcf` containing the variations detected between the reference genome and the sample KMC DB.

## 4. Find IBS windows
To identify IBS (Identity by State) windows, use the `findIBS` command:

    $ kcftools findIBS \
        --input sample.kcf \
        --output sample.ibs \
        --summary

This command will generate an IBS summary file named `sample.ibs.summary.tsv` containing the IBS windows identified in the KCF file with additional summary.

!!! note
    A workaround [bash script](https://github.com/sivasubramanics/kcftools/blob/main/utils/run_kcftools.sh) to run the complete pipeline for multiple samples is provided along with the KCFTOOLS source code. You can find it in the `utils` directory. This script automates the process of creating KCF files for multiple samples and finding IBS windows. Details on how to use this script can be found in the [Run Pipeline](usage/pipeline.md) documentation.


