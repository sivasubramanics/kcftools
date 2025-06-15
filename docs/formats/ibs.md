# IBS Summary File

This document describes the structure of the **IBS Summary File**, a key output used in comparative genomics analyses to record *Identical By State* (IBS) segments across genomic intervals for a given sample. Along with this summary file, the `findIBS` command can also produce a BED file format output, which is useful for visualization in genome browsers.

---

## File Format

The IBS summary file is a **tab-delimited text file** where each row represents a genomic region (blocks) for a given sample. Each column is described in detail below.

### Column Definitions

| Column         | Data Type | Description |
|----------------|------------|-------------|
| `Block`        | `int`      | Unique identifier for the genomic block (can be sequential or tool-specific). |
| `Sample`       | `string`   | Sample identifier. It is recommended to use standardized identifiers (e.g., `SAM001`, `SAM002`) instead of internal IDs like `LK067`. |
| `Chromosome`   | `string`   | Name of the chromosome or scaffold. Standard naming conventions such as `chr1`, `chr2`, ..., `chrUn` (for unplaced scaffolds) should be used. |
| `Start`        | `int`      | Start position (0-based, inclusive) of the block in base pairs. |
| `End`          | `int`      | End position (0-based, inclusive) of the block in base pairs. |
| `Length`       | `int`      | Length of the block in base pairs. Calculated as `End - Start + 1`. |
| `TotalBlocks`  | `int`      | Total number of sub-units (or smaller analysis windows) within this block. |
| `IBSBlocks`    | `int`      | Number of sub-units within the block that were classified as IBS. |
| `IBSProportion`| `float`    | Fraction of sub-blocks that are IBS. Computed as `IBSBlocks / TotalBlocks`. Ranges from `0.0` to `1.0`. |
| `MeanScore`    | `float`    | Mean similarity score across the block (can reflect allele sharing, IBS strength, or other scoring mechanisms depending on the tool used). |

---

## Example

```text
Block	Sample	Chromosome	Start	End	Length	TotalBlocks	IBSBlocks	IBSProportion	MeanScore
3495	SAM067	chrUn_01	0	15194	15194	4	4	1.00	0.03
3475	SAM067	chrUn_02	0	59659	59659	12	9	0.75	74.59
3148	SAM067	chr5	1421134	1426134	5000	1	1	1.00	92.11
```

---

!!! note
    - A **high IBSProportion** (`1.00`) indicates that the entire region shows complete similarity across all tested windows.
    - A **low MeanScore** may indicate sparse or weak similarity, even if all windows are IBS.
    - A **block with IBSProportion < 1.0** can be indicative of genomic heterogeneity or recombination signals in that region.
    - **Short blocks with high IBS** might be less informative than longer blocks with consistent similarity across a large region.

---