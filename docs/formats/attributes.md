# Attribute Files

Attribute files are just a simple tab-separated values (TSV) file that contains attributes for each atrributes (eg: observed *k*-mer, variations, score, etc). Each line corresponds to a window, and each columns corresponds to a sample. The first column is the window ID, and the subsequent columns are sample-specific attributes.

| Field         | Description                                           |
|---------------|-------------------------------------------------------|
| `window_id` | Unique identifier for the window                      |
| `SAM001` | Sample-specific attribute value (e.g., observed *k*-mers, variations, etc.) |
| `SAM002` | Sample-specific attribute value (e.g., observed *k*-mers, variations, etc.) |
| `SAM003` | Sample-specific attribute value (e.g., observed *k*-mers, variations, etc.) |

There are numerous attribute files that will be generated by the `getAttributes` command, each corresponding to a specific attribute. The attributes are defined in the KCF file format and can include:

| Attribute File Exetension | Description                                           |
|---------------------------|-------------------------------------------------------|
| `.obs.tsv`                | Observed *k*-mers for each sample in each window      |
| `.var.tsv`                | Variations detected in each window for each sample    |
| `.score.tsv`              | Score for each sample in each window                  |
| `.inDist.tsv`             | Inner distance for each sample in each window         |
| `.tailDist.tsv`           | Tail distance for each sample in each window          |
| `.totalkmers.tsv`         | Total *k*-mers in each window                         |

#### Example Attribute File

| window_id     | SAM001 | SAM002 | SAM003 | SAM004 |
|---------------|--------|--------|--------|--------|
| chr1_0        | 81.50  | 81.42  | 81.50  | 81.49  |
| chr1_4969     | 67.04  | 67.09  | 66.75  | 66.99  |
| chr1_9938     | 73.60  | 73.62  | 74.71  | 89.61  |
| chr1_14907    | 62.50  | 62.50  | 62.50  | 87.11  |
| chr1_19876    | 76.26  | 76.26  | 65.64  | 76.28  |
| chr1_24845    | 45.46  | 45.52  | 45.57  | 45.90  |
| chr1_29814    | 90.06  | 89.84  | 76.59  | 89.22  |
| chr1_34783    | 91.86  | 92.76  | 75.22  | 92.57  |
| chr1_39752    | 93.44  | 93.41  | 91.71  | 93.36  |



