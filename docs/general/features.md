# Features

KCFTOOLS provides a comprehensive set of features for genomic analysis, leveraging *k*-mer counting techniques to facilitate various tasks. Below is an overview of the key features:

#### [**getVariations**](usage/getVariations.md)
- Identifies variations between a reference genome and a query sample by counting *k*-mer presence/absence.
- Supports different windowing strategies (fixed-length, gene models, or transcript features).
- Outputs a KCF file containing the variations with detailed identity scores. 
- Allows memory-efficient processing with options for multi-threading.
   
#### [**findIBS**](usage/findIBS.md)
- Identifies **Identity-by-State (IBS)** windows or **Variable Regions (VR)** in a KCF file.
- Calculates identity scores for each window and filters them based on user-defined thresholds.
- Merges consequent windows to create larger regions of identical/variable regions.
- Outputs a summary file with detailed information about the IBS/VR windows.

#### [**cohort**](usage/cohort.md)
- Combines multiple KCF files (sample-wise) into a single cohort file.
- Facilitates cohort-level analyses by aggregating variations across multiple samples.
- Maintains the number of windows and reference sequence information for each sample.

#### [**kcfToMatrix**](usage/kcfToMatrix.md)
- Converts a KCF file into a genotype matrix format.
- Outputs a matrix where rows represent samples and columns represent allele codes (0, 1, 2).
- Outputs a map file that maps the window IDs to their respective positions in the reference genome.
- The output files are suitable for further analyses in tools like [PLINK](https://www.cog-genomics.org/plink/), [Tassel](https://www.maizegenetics.net/tassel) or [GAPIT](https://zzlab.net/GAPIT/)
- Supports options for filtering based on allele frequency and minimum missing data.

#### [**kcf2tsv**](usage/kcf2tsv.md)
- Converts a KCF file into a tab-separated values (TSV) format.
- Outputs a file with detailed information about each window, including start and end positions, total *k*-mers, observed *k*-mers, variations and kmer_distance (identical to IBSpy output).

#### [**getAttributes**](usage/getAttributes.md)
- Extracts attributes from a KCF file and outputs them in a tab-separated format.
- Helps in retrieving specific information about the variations, such as *k*-mer counts and identity scores.
- This data could be used effectively for plotting and visualization of variation data.

#### [**splitKCF**](usage/splitKCF.md)
- Splits a KCF file into multiple smaller files based on chromosome or contig.
- Facilitates parallel processing of large KCF files by breaking them down into manageable chunks.
- Maintains the integrity of the KCF file format while allowing for efficient data handling.

#### [**increaseWindow**](usage/increaseWindow.md)
- Increases the window size of a KCF file by a specified factor.
- Useful for adjusting the resolution of the variations detected in the KCF file.
- Maintains the original KCF file structure while expanding the analysis window.
- Supports options for adjusting the effective length and identity score calculations based on the new window size.
- Allows for re-evaluation of variations with the new window size, potentially avoiding the need to re-run the `getVariations` command.
- Outputs a new KCF file with the adjusted window size and recalculated identity scores.
