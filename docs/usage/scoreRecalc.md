# kcftools scoreRecalc

The `kcftools scoreRecalc` command is used to recalculate identity scores in a KCF file based on updated weights. This allows users to adjust scoring criteria dynamically without the need to re-run the entire variation detection process.

---
## Usage

    $ kcftools scoreRecalc -i=<inFile> -o=<outFile>
                            [--wi=<innerDistanceWeight>]
                            [--wr=<kmerRatioWeight>] [--wt=<tailDistanceWeight>]

---

## Description
The `scoreRecalc` command processes an existing KCF file and recalculates the identity scores for each window using new weight parameters. This is particularly useful when users want to refine their scoring criteria based on new insights or requirements without having to repeat the computationally intensive steps of variation detection.

## Options
| Option                                | Description                                                                 | Default    | Required |
|---------------------------------------|-----------------------------------------------------------------------------|------------|----------|
| `-i`, `--input=<inFile>`              | Input `.kcf` file to recalculate scores                                  | _N/A_      | Yes      |
| `-o`, `--output=<outFile>`          | Output `.kcf` file with recalculated scores                              | _N/A_      | Yes      |
| `--wi=<innerDistanceWeight>`            | Weight for inner distance in score calculation                        | `0.5`     | No       |
| `--wr=<kmerRatioWeight>`            | Weight for k-mer ratio in score calculation                        | `0.3`     | No       |
| `--wt=<tailDistanceWeight>`            | Weight for tail distance in score calculation                        | `0.2`     | No       |
| `--help`                              | Show help message and exit                                                | _N/A_      | No       |

---
## Output
This command generates a new KCF file with recalculated identity scores based on the specified weights.

**Example output file:**

    recalculated_scores.kcf
---

## Example
**Recalculate scores in a KCF file with custom weights:**

    $ kcftools scoreRecalc -i original.kcf -o recalculated_scores.kcf --wi 0.6 --wr 0.2 --wt 0.2

This command:
- Uses an inner distance weight of `0.6`
- Uses a k-mer ratio weight of `0.2`
- Uses a tail distance weight of `0.2`
- Generates `recalculated_scores.kcf` with updated identity scores based on the new weights

---

## Help
For more information on the `scoreRecalc` command and its options, run:
    
    $ kcftools scoreRecalc --help
