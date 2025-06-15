# KCF File Format

The **Kmer Count Format (`.kcf`)** file summarizes the variation profile of a query relative to a reference genome using *k*-mer presence/absence matrices.

---

## KCF File Header Description

The `.kcf` file begins with a series of metadata headers describing the format version, reference genome, software used, and analysis parameters.

### Example Header Lines

```
##format=KCF0.1
##date=2024-12-05
##source=kcftools
##reference=lsatv11.chr3.fasta
##contig=<ID=chr3,length=324658466>
##INFO=<ID=IS,Type=Float,Description="Minimum score for the window">
##INFO=<ID=XS,Type=Float,Description="Maximum score for the window">
##INFO=<ID=MS,Type=Float,Description="Mean score for the window">
##INFO=<ID=IO,Type=Integer,Description="Minimum observed kmers in the window">
##INFO=<ID=XO,Type=Integer,Description="Maximum observed kmers in the window">
##INFO=<ID=MO,Type=Integer,Description="Mean observed kmers in the window">
##INFO=<ID=IV,Type=Integer,Description="Minimum variations in the window">
##INFO=<ID=XV,Type=Integer,Description="Maximum variations in the window">
##INFO=<ID=MV,Type=Integer,Description="Mean variations in the window">
##FORMAT=<ID=IB,Type=Integer,Description="IBS number">
##FORMAT=<ID=VA,Type=Integer,Description="Variations">
##FORMAT=<ID=OB,Type=Integer,Description="Observed kmers">
##FORMAT=<ID=ID,Type=Integer,Description="Inner distance">
##FORMAT=<ID=LD,Type=Integer,Description="Left tail distance">
##FORMAT=<ID=RD,Type=Integer,Description="Right tail distance">
##FORMAT=<ID=SC,Type=Float,Description="Score">
##PARAM=<ID=window,value=50000>
##PARAM=<ID=kmer,value=31>
##PARAM=<ID=IBS,value=false>
##PARAM=<ID=nwindow,value=6498>
##CMD=kcftools-0.0.1-SNAPSHOT.jar getVariations -k lsal.chr3 -o lsal.kcftools.kcf -r lsatv11.chr3.fasta -s lsal -t 24 -w 50000
```

### Header Line Details

| Header Key     | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `##format`     | Version of the KCF format specification.                                    |
| `##date`       | Date of generation.                                                         |
| `##source`     | Software/tool that created the file.                                        |
| `##reference`  | Reference genome FASTA file used in analysis.                               |
| `##contig`     | Reference contig ID and its length.                                         |
| `##INFO`       | Window-level summary statistics (score, kmers, variations).                 |
| `##FORMAT`     | Sample-level metrics (e.g., IB, VA, OB, etc.)                              |
| `##PARAM`      | Runtime parameters used in the tool.                                        |
| `##CMD`        | Full command used to generate the KCF file (for reproducibility).           |

---

## KCF File Data Columns

Each row of the file (after headers) represents a **non-overlapping genomic window** with its variation summary.

### Main Columns

| Column         | Description                                                                 |
|----------------|-----------------------------------------------------------------------------|
| `CHROM`        | Reference chromosome or contig name.                                        |
| `START`        | Start coordinate of the window (0-based).                                   |
| `END`          | End coordinate of the window (0-based, exclusive).                          |
| `TOTAL_KMERS`  | Total reference *k*-mers in the window.                                     |
| `INFO`         | Semicolon-separated list of window-level stats (IS, XS, MS, etc.).          |
| `FORMAT`       | Format string that defines the sample field structure (e.g., IB:VA:OB:...). |
| `<Sample>`     | Sample-specific colon-separated values based on the `FORMAT`.               |

---

### Format Field Attributes

Sample fields are defined based on the `FORMAT` column and provide detailed metrics per window.

| Field | Description                                                                |
|-------|----------------------------------------------------------------------------|
| `IB`  | Identity-by-state number (shared *k*-mers)                                 |
| `VA`  | Number of variations detected in the window                               |
| `OB`  | Number of observed *k*-mers (sample-specific presence)                     |
| `ID`  | Inner distance (gap within the window between *k*-mer hits)                |
| `LD`  | Left tail distance (gap at the beginning of the window)                    |
| `RD`  | Right tail distance (gap at the end of the window)                         |
| `SC`  | Identity score (e.g., percentage similarity between sample and reference)  |

---

!!! note
    - All coordinates follow the **0-based** BED-style convention.
    - The `FORMAT` field is essential for decoding sample data.
    - INFO fields are useful for filtering or plotting summaries at a per-window level.

---

## Example Data Row

```
chr3    0       50000   1012    IS=0.91;XS=1.00;MS=0.95;IO=980;XO=1000;MO=990;IV=0;XV=3;MV=1     IB:VA:OB:ID:LD:RD:SC     989:0:990:0:0:0:99.95
```

This row shows:
- A 50kb window on chromosome 3
- 1012 *k*-mers total
- High identity score (99.95)
- 0 variations detected
- Sample-specific field values are provided in the last column

