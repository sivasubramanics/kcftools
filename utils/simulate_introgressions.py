#!/usr/bin/env python3
# Script: simulate_introgressions.py
# Description: Simulate introgressions from a donor genome into a recurrent genome.
# Note: Introgression size and numbers are adjustable.
# Author: c.s.sivsubramani@gmail.com
# Date: 2025-09-23

"""
This script takes two fasta files and introduces introgressions from a donor genome into a recurrent genome.
It also returns the locations of the introgressions, updates GTF files accordingly, and generates a visualization plot.
"""

import argparse
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class FASTA:
    def __init__(self, name: str, sequence: str, description: str = None):
        self.name = name
        self.sequence = sequence
        self.description = description

    def __str__(self):
        return self.fold(60)

    def fold(self, width: int) -> str:
        out_str = f">{self.name}"
        if self.description:
            out_str += f" {self.description}"
        out_str += "\n"
        out_str += '\n'.join([self.sequence[i:i + width] for i in range(0, len(self.sequence), width)])
        return out_str

    def __len__(self):
        return len(self.sequence)

def parse_fasta(file):
    name = None
    description = None
    sequence = ""
    with open(file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield FASTA(name, sequence, description)
                name = line[1:].split()[0]
                description = line[len(name) + 2:]
                sequence = ""
            else:
                sequence += line
        yield FASTA(name, sequence, description)

def introduce_introgressions(donor_seq, recurrent_seq, num_introgressions=2,
                             donor_size_range=(1_000_000, 3_000_000),
                             recipient_size_range=(1_000_000, 3_000_000),
                             max_attempts=1000):
    donor_len = len(donor_seq.sequence)
    recurrent_len = len(recurrent_seq.sequence)

    introgressions_applied = []
    introgressed_seq = recurrent_seq.sequence
    offset = 0

    # --- Step 1: choose donor intervals (non-overlapping, sorted) ---
    donor_intervals = []
    for _ in range(num_introgressions):
        attempt = 0
        while attempt < max_attempts:
            attempt += 1
            donor_size = random.randint(*donor_size_range)
            d_start = random.randint(0, donor_len - donor_size)
            d_end = d_start + donor_size
            if all(d_end <= u[0] or d_start >= u[1] for u in donor_intervals):
                donor_intervals.append((d_start, d_end))
                break
        else:
            raise RuntimeError("Could not place non-overlapping donor interval")
    donor_intervals.sort()

    # --- Step 2: choose recurrent intervals (non-overlapping, sorted) ---
    recurrent_intervals = []
    for _ in range(num_introgressions):
        attempt = 0
        while attempt < max_attempts:
            attempt += 1
            recipient_size = random.randint(*recipient_size_range)
            r_start = random.randint(0, recurrent_len - recipient_size)
            r_end = r_start + recipient_size
            if all(r_end <= u[0] or r_start >= u[1] for u in recurrent_intervals):
                recurrent_intervals.append((r_start, r_end))
                break
        else:
            raise RuntimeError("Could not place non-overlapping recurrent interval")
    recurrent_intervals.sort()

    # --- Step 3: apply introgressions in order ---
    for (d_start, d_end), (r_start, r_end) in zip(donor_intervals, recurrent_intervals):
        d_frag = donor_seq.sequence[d_start:d_end]
        adj_r_start = r_start + offset
        adj_r_end = r_end + offset

        introgressed_seq = (
            introgressed_seq[:adj_r_start] +
            d_frag +
            introgressed_seq[adj_r_end:]
        )

        introgressions_applied.append({
            "donor_start": d_start,
            "donor_end": d_end,
            "recurrent_start": r_start,
            "recurrent_end": r_end,
            "new_start": adj_r_start,
            "new_end": adj_r_start + len(d_frag),
            "donor_length": len(d_frag),
            "length_diff": len(d_frag) - (r_end - r_start)
        })

        offset += len(d_frag) - (r_end - r_start)

    return introgressed_seq, introgressions_applied



def update_gtf(gtf_file, chrom_name, new_chrom, introgressions):
    updated_lines = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if fields[0] != chrom_name:
                continue
            start = int(fields[3])
            end = int(fields[4])
            for reg in introgressions:
                if start >= reg['donor_start'] and end <= reg['donor_end']:
                    new_start = reg.get('new_start', reg['new_end'] - (reg['donor_end'] - reg['donor_start']))
                    shift = new_start - reg['donor_start']
                    fields[0] = new_chrom
                    fields[3] = str(start + shift)
                    fields[4] = str(end + shift)
                    updated_lines.append("\t".join(fields))
                    break
    return updated_lines


def extract_non_introgressed_gtf(gtf_file, chrom_name, new_chrom, introgressions):
    updated_lines = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if fields[0] != chrom_name:
                continue
            start = int(fields[3])
            end = int(fields[4])
            shift = 0
            include = True
            for reg in introgressions:
                if start < reg['recurrent_end'] and end > reg['recurrent_start']:
                    include = False
                    break
                elif start >= reg['recurrent_end']:
                    shift += reg['length_diff']
            if include:
                fields[0] = new_chrom
                fields[3] = str(start + shift)
                fields[4] = str(end + shift)
                updated_lines.append("\t".join(fields))
    return updated_lines

def draw_introgression_plot(donor_len, recurrent_len, introgressed_len, introgressions, output_prefix):
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.set_xlim(0, max(donor_len, recurrent_len, introgressed_len))
    ax.set_ylim(0, 6)
    ax.axis("off")

    def draw_band(y, length, label, color):
        rect = patches.Rectangle((0, y), length, 0.4, linewidth=1, edgecolor='black', facecolor=color)
        ax.add_patch(rect)
        ax.text(length + 5_000_000, y + 0.2, label, va='center', fontsize=10)

    donor_color = 'skyblue'
    recurrent_color = 'wheat'

    draw_band(5, donor_len, "Donor Genome", donor_color)
    draw_band(3, introgressed_len, "Introgressed Genome", recurrent_color)
    draw_band(1, recurrent_len, "Recurrent Genome", recurrent_color)

    for intro in introgressions:
        donor_len = intro.get("donor_length", intro["donor_end"] - intro["donor_start"])
        new_start = intro.get("new_start", intro["new_end"] - donor_len)
        donor_mid = (intro["donor_start"] + intro["donor_end"]) / 2
        new_mid = (new_start + intro["new_end"]) / 2
        intro_length = donor_len
        ax.plot([donor_mid, new_mid], [5, 3.4], linestyle='--', alpha=0.8)
        intro_rect = patches.Rectangle((new_start, 3), intro_length, 0.4,
                                        linewidth=1, edgecolor='black', facecolor=donor_color)
        ax.add_patch(intro_rect)
        donor_rect = patches.Rectangle((intro["donor_start"], 5),
                                        intro["donor_end"] - intro["donor_start"], 0.4,
                                        linewidth=2, edgecolor='black', facecolor=donor_color)
        ax.add_patch(donor_rect)
        removed_rect = patches.Rectangle((intro["recurrent_start"], 1),
                                            intro["recurrent_end"] - intro["recurrent_start"], 0.4,
                                            linewidth=1, edgecolor='red', facecolor='none', linestyle='--')
        ax.add_patch(removed_rect)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_introgressed.svg")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Introduce introgressions and output sequence, annotation and plot')
    parser.add_argument('-d', '--donor', required=True)
    parser.add_argument('-e', '--donor_gtf', required=True)
    parser.add_argument('-r', '--recurrent', required=True)
    parser.add_argument('-s', '--recurrent_gtf', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-n', '--num_introgressions', type=int, default=2)
    args = parser.parse_args()

    donor_seq = next(parse_fasta(args.donor))
    recurrent_seq = next(parse_fasta(args.recurrent))

    new_seq_str, introgressions = introduce_introgressions(donor_seq, recurrent_seq, args.num_introgressions)
    introgressed_fasta = FASTA("intr", new_seq_str, "Introgressed chromosome")

    with open(f"{args.output}_introgressed.fasta", "w") as f:
        f.write(str(introgressed_fasta) + "\n")

    with open(f"{args.output}_introgressed.tsv", "w") as f:
        f.write("Donor_Chr\tDonor_Start\tDonor_End\tRecurrent_Chr\tRecurrent_Start\tRecurrent_End\tLength\n")
        for reg in introgressions:
            f.write(f"{donor_seq.name}\t{reg['donor_start']}\t{reg['donor_end']}\t"
                    f"{recurrent_seq.name}\t{reg['recurrent_start']}\t{reg['recurrent_end']}\t"
                    f"{reg['donor_end'] - reg['donor_start']}\n")

    gtf_lines = []
    gtf_lines.extend(extract_non_introgressed_gtf(args.recurrent_gtf, recurrent_seq.name, "intr", introgressions))
    gtf_lines.extend(update_gtf(args.donor_gtf, donor_seq.name, "intr", introgressions))

    with open(f"{args.output}_introgressed.gtf", "w") as f:
        for line in gtf_lines:
            f.write(line + "\n")

    draw_introgression_plot(len(donor_seq), len(recurrent_seq), len(new_seq_str), introgressions, args.output)
    print("All files written successfully.")

if __name__ == "__main__":
    main()
