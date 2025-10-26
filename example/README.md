# Detecting Simulated Introgressions using **KCFtools**

This example workflow demonstrates how to use **KCFtools** to detect introgressed regions from simulated genomic data.  

---

## Example Input Files

This directory contains three FASTA files and two additional files:

| File           | Description                                                         |
|----------------|---------------------------------------------------------------------|
| `don.fasta.gz` | Donor genome (introgression source)                                 |
| `rec.fasta.gz` | Recipient genome (recurrent parent)                                 |
| `itr.fasta.gz` | Simulated introgressed genome (donor segments introgressed into recipient) |
| `itr.tsv`     | True introgressed segments information                              |
| `itr.svg`     | Visualization of true introgressed segments                        |

We used `utils/simulate_introgressions.py` script to simulate the introgressions.

`itr.svg` and `itr.tsv` shows the true introgressed segments in `itr.fasta.gz`.

---

## Prerequisites

Ensure you have installed `kcftools` and `KMC` and they are accessible.

Note: The R libraries are required for plotIBS.R is not included in the conda environment. You may need to install them separately.

#### Uncompress the input FASTA files:

```bash
gzip -d -k don.fasta.gz rec.fasta.gz itr.fasta.gz
```
This will create `don.fasta`, `rec.fasta`, and `itr.fasta` files.

---

## Step-by-Step Workflow

### **1. Count k-mers from the introgressed genome**

```bash
kmc -k31 -m4 -t4 -ci0 -p9 -fm itr.fasta itr.k31 .
```
This command counts 31-mers from the `itr.fasta` file and stores the results in `itr.k31` database.

parameters:
- `-k31`: k-mer size of 31
- `-m4`: use up to 4GB of RAM
- `-t4`: use 4 threads
- `-ci0`: count all k-mers (no minimum count)
- `-p9`: use signature length of 9
- `-fm`: input is in FASTA format
- `itr.fasta`: input FASTA file
- `itr.k31`: output KMC database name
- `.`: current directory for temporary files

---

### **2. Generate KCF file comparing introgressed genome to donor genome**

```bash
kcftools getVariations -r don.fasta -k itr.k31 -o don.itr.k31.w5k.kcf -s itr -f window -w 5000
```
This command generates a KCF file `don.itr.k31.w5k.kcf` comparing the introgressed genome to the donor genome using 5kb windows.

parameters:
- `getVariations`: KCFtools command to get variations
- `-r don.fasta`: reference genome (donor)
- `-k itr.k31`: KMC database of introgressed genome
- `-o don.itr.k31.w5k.kcf`: output KCF file
- `-s itr`: sample name for introgressed genome
- `-f window`: use window-based analysis
- `-w 5000`: window size of 5000 bp

### **3. Identify introgressed regions using findIBS**

```bash
kcftools findIBS -i don.itr.k31.w5k.kcf -o don.itr.k31.w5k.ibs --min 8 --score 90 --summary
```

This command identifies introgressed regions by finding IBS segments in the KCF file.

parameters:
- `findIBS`: KCFtools command to find IBS segments
- `-i don.itr.k31.w5k.kcf`: input KCF file
- `-o don.itr.k31.w5k.ibs`: output IBS file prefix
- `--min 8`: minimum number of markers in an IBS segment
- `--score 90`: minimum IBS score threshold
- `--summary`: generate a summary of IBS segments

### **4. Plot the detected introgressed regions**
```bash
# Generate chromosome info file
awk '{i++; s += $2; printf "%s\t%s\t%d\t%.0f\n", $1, $2, i, s }' don.fasta.faidx > don.chrom_info.tsv
# Plot IBS segments
Rscript utils/plotIBS.R -c don.chrom_info.tsv -i don.itr.k31.w5k.ibs.summary.tsv -o don.itr.k31.w5k.ibs.pdf
```
This command generates a plot of the detected introgressed regions using the IBS summary file. Make sure the plotIBS.R script path is correct.

parameters:
- `-c don.chrom_info.tsv`: chromosome information file
- `-i don.itr.k31.w5k.ibs.summary.tsv`: IBS summary file
- `-o don.itr.k31.w5k.ibs.pdf`: output PDF file

---

### **5. Generate KCF file comparing introgressed genome to recipient genome**

```bash
kcftools getVariations -r rec.fasta -k itr.k31 -o rec.itr.k31.w5k.kcf -s itr -f window -w 5000
```

This command generates a KCF file `rec.itr.k31.w5k.kcf` comparing the introgressed genome to the recipient genome using 5kb windows.
parameters:
- `-r rec.fasta`: reference genome (recipient)
- `rec.itr.k31.w5k.kcf`: output KCF file
- other parameters are the same as in step 2

---

### **6. Identify Variable regions using findIBS**

```bash
kcftools findIBS -i rec.itr.k31.w5k.kcf -o rec.itr.k31.w5k.ibs --min 8 --score 90 --summary --var
```
This command identifies variable regions by finding IBS segments in the KCF file.
parameters:
- `-i rec.itr.k31.w5k.kcf`: input KCF file
- `-o rec.itr.k31.w5k.ibs`: output IBS file prefix
- `--var`: identify variable regions
- other parameters are the same as in step 3

### **7. Plot the detected variable regions**
```bash
# Generate chromosome info file
awk '{i++; s += $2; printf "%s\t%s\t%d\t%.0f\n", $1, $2, i, s }' rec.fasta.faidx > rec.chrom_info.tsv

# Plot variable segments
Rscript utils/plotIBS.R -c rec.chrom_info.tsv -i rec.itr.k31.w5k.ibs.summary.tsv -o rec.itr.k31.w5k.ibs.pdf --var
```
This command generates a plot of the detected variable regions using the IBS summary file. Make sure the plotIBS.R script path is correct.

parameters:
- `-c rec.chrom_info.tsv`: chromosome information file
- `-i rec.itr.k31.w5k.ibs.summary.tsv`: IBS summary
- `-o rec.itr.k31.w5k.ibs.pdf`: output PDF file
- `--var`: indicate variable regions
- other parameters are the same as in step 4

---
## Results Interpretation
- The PDF files `don.itr.k31.w5k.ibs.pdf` and `rec.itr.k31.w5k.ibs.pdf` visualize the detected introgressed and variable regions, respectively.
- Compare these results with the true introgressed segments in `itr.svg` and `itr.tsv` to evaluate the accuracy of detection.
- Adjust parameters in `getVariations` and `findIBS` as needed to optimize detection based on your specific dataset and requirements.
