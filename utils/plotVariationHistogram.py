#!/usr/bin/env python3

"""
This script is used to plot the histogram of the variation counts per window
"""

import argparse
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter

def get_variation_counts(kcf_file, sample=None):
    with open(kcf_file) as f:
        variation_counts = []
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                samples = line.strip().split('\t')[6:]
                continue
            fields = line.strip().split('\t')
            if sample is not None and sample not in samples:
                logging.error(f'Sample {sample} not found in KCF file')
                return
            if sample is not None:
                sample_index = samples.index(sample)
                # get variation count from sample column GT:VA:OB:SC as int list
                variation_count = [int(field.split(':')[1]) for field in fields[6+sample_index].split(',')]
            else:
                # get variation count from all samples as int list
                variation_count = [int(field.split(':')[1]) for field in fields[6].split(',')]
            variation_counts.extend(variation_count)
    return variation_counts

# def plot_histogram(df, score, fig_size):
#     font_s = 15
#     font_w = "bold"
#     fig = plt.figure(figsize=fig_size, facecolor='w')
#     #     plt.ylim(0,20000)
#     plt.xscale('symlog', base=2)
#     bin_range = range(min(df[score]), max(df[score]))
#     print(bin_range)
#     sns.histplot(data=df, x=score, kde=True, bins=bin_range)
#     plt.xticks([-0.55, 0, 1, 2, 3, 5, 10, 20, 40, 100, 300, 1000, 2000], fontsize=12)
#     plt.yticks(fontsize=font_s)
#     fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
#     plt.title(f'reference: {reference}, query: {query}, window: {window}')
#     plt.xlabel(score, fontsize=font_s)
#     plt.ylabel('Count', fontsize=font_s, fontweight=font_w)
#     sns.despine(offset=5)


def plot_histogram(df, score, fig_size, output):
    font_s = 15
    font_w = "bold"
    fig = plt.figure(figsize=fig_size, facecolor='w')
    plt.xscale('symlog', base=2)

    bin_range = range(min(df[score]), max(df[score]))
    print(bin_range)

    # Use the numpy array in histplot
    sns.histplot(data=df, x=score, kde=True, bins=bin_range)

    plt.xticks([-0.55, 0, 1, 2, 3, 5, 10, 20, 40, 100, 300, 1000, 2000], fontsize=12)
    plt.yticks(fontsize=font_s)
    fig.axes[0].xaxis.set_major_formatter(ScalarFormatter())
    plt.title(f'reference: {reference}, query: {query}, window: {window}')
    plt.xlabel(score, fontsize=font_s)
    plt.ylabel('Count', fontsize=font_s, fontweight=font_w)
    sns.despine(offset=5)
    plt.savefig(output, dpi=300)



def main():
    parser = argparse.ArgumentParser(description='Plot histogram of variation counts per window from KCF file')
    parser.add_argument('-i', '--input', help='Input KCF file', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

    # variation_counts = get_variation_counts(args.input)
    variation_counts = pd.read_csv(args.input, sep='\t')
    print(variation_counts.head())
    if variation_counts is None:
        logging.error('Error getting variation counts')
        return

    plot_histogram(variation_counts, "LK067", (25, 4), args.output)
    logging.info(f'Output saved to {args.output}')


if __name__ == '__main__':
    main()