# kcftools
This repository holds scripts that can be used to detect introgressions form the resequencing data

## detect variations
### IBSpy_call_variations.sh
shell script can be used to run the fastq->kmer_variants using the KMC and IBSpy pipeline. Before running the pipeline, make sure the dependencies are installed and added to the env path

`Usage: IBSpy_call_variation.sh [-h] [-f fq.list] [-F fq_list] [-r ref.fa] [-o output] [-p threads] [-k kmer_size] [-w window_size]
Options:
  -f STR    list of fastq files
  -F STR    format (fq/fq_list/fa/fa_list/)
  -r STR    reference genome
  -o STR    output file
  -p INT    number of threads
  -k INT    kmer size
  -w INT    window size
  -h        help
  `

#### dependencies
[KMC](https://github.com/refresh-bio/KMC)  
[kmerGWAS](https://github.com/voichek/kmersGWAS)  
[IBSpy_cpp](https://github.com/Uauy-Lab/IBSpy)  

## kcftools 

`usage: kcftools.py [-h] {tsv2kcf,cohort,increase_window,find_IBS,extract,kcf2bedgraph,kcf2matrix,split_kcf,concat} ...

positional arguments:
  {tsv2kcf,cohort,increase_window,find_IBS,extract,kcf2bedgraph,kcf2matrix,split_kcf,concat}
    tsv2kcf             Step 1: Convert IBSpy tsv to kcf
    cohort              Step 2.a: Combine kcf files
    increase_window     Step 2.b: Increase window size
    find_IBS            Step 3: Find the IBS regions in kcf file
    extract             Step 4: Extract windows from kcf file
    kcf2bedgraph        Step 5: Convert kcf file to bedgraph files
    kcf2matrix          Step 6: Convert kcf file to genotype matrix
    split_kcf           Misc: Split kcf file by chromosome
    concat              Misc: Concatenate kcf files (samples should be identical)

optional arguments:
  -h, --help            show this help message and exit
  `

### tsv2kcf 
inputs: single TSV file  
convert IBSpy TSV file to KCF format. This can be used for only one TSV file at a time

### cohort 
inputs: multiple KCF files  
combine multiple KCF files and output a multi sample KCF file

### increase_window 
inputs: single/multisample KCF, window length, kmersize  
To change the windows from the KCF file. This function will add the consecutive variation counts for the merged windows

### find_IBS
inputs: multisample KCF file  
Based on the score (OBS-VAR/TOTAL), this script will predict if the given window is an IBS regions. Also this will check for a consecutive windows to have IBS and call it an possible introgression

### extract
inputs: find_IBS output KCF file  
It will extract the possible introgression blocks form the find_IBS output

### kcf2bedgraph
inputs: find_IBS output KCF file  
Will convert find_IBS output KCF to bedgraph, which can be used in any GenomeBrowsers view the introgressions (or IBS)

### kcf2matrix 
inputs: multisample KCF  
will convert the KCF format to the numberical genotype format, so that can be used for GWAS or QTL analysis

### split_kcf 
inputs: single/multisample KCF  
will split the KCF files to chromosome wise

### concat 
inputs: single/multi sample KCFs  
will combine the KCF for different chromosomes to one KCF file. The sample count and order MUST be identical between the input files.


please reach siva.selvanayagam[at]wur.nl for bugfixes and features
