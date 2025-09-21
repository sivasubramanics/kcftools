# Genotype Data Table

The KCF Genotype Table file format is designed to represent the genotype data extracted from KCF files in a tabular format suitable for further analysis in tools like Tassel, or GAPIT.

The Genotype data consists of three files:

1. **Genotype Table**: A table where columns represent samples and rows represent allele codes (0, 1, 2, -1).
2. **Contigs Map File**: A file that maps contig names to their respective numeric id.

## Genotype Table Format
The table file is a tab-separated values (TSV) file that contains metadata about the regions (windows) and alleles. Each line corresponds to a window and includes the following fields:

| Field   | Description                                           |
|---------|-------------------------------------------------------|
| `ID`    | Unique identifier for the window                      |
| `CHR`   | Chromosome or contig name                             |
| `START` | Start position of the window (0-based)                |
| `END`   | End position of the window (1-based, exclusive)       |


#### Example Genotype Table File

| ID | CHR |START|	END|	SAM01|	SAM02|	SAM03|	SAM04|	SAM05|	SAM06|	SAM07|	SAM08|	SAM09|	SAM10|	SAM11|
|------|------|------|------|------|------|------|------|------|------|------|------|------|------|------|
|1|	1|	0|	1000|	0|	0|	0|	0|	0|	0|	0|	0|	0|	0|	0|
|2|	1|	970|	1970|	0|	0|	0|	0|	0|	0|	0|	0|	0|	0|	0|
|3|	1|	1940|	2940|	0|	0|	0|	0|	0|	0|	0|	0|	0|	0|	0|
|4|	1|	2910|	3910|	-1|	-1|	-1|	-1|	-1|	-1|	-1|	-1|	-1|	-1|	-1|
|5|	1|	3880|	4880|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|
|6|	1|	4850|	5850|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|
|7|	1|	5820|	6820|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|
|8|	1|	6790|	7790|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|
|9|	1|	7760|	8760|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|	2|


## Contigs Map File Format
The contigs map file is a tab-separated values (TSV) file that maps contig names to their respective numeric id.

| contig_name | contig_id |
|--------------|-----------|
| chr1         | 1         |
| chr2         | 2         |

