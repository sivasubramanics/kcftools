#!/usr/bin/env python3

"""
functions for writing parsers files
"""
import csv
import logging
import sys

from kcftools.utils.helper_functions import monomorphic


def write_kcf(out_kcf, windows, samples, misc_lines=None):
    """
    Write kcf file
    """
    if type(samples) == str:
        samples = [samples]
    logging.info(f'Writing {out_kcf}')
    with open(out_kcf, 'w') as o:
        samples_line = '\t'.join(samples)
        if misc_lines:
            for misc_line in misc_lines:
                o.write(f'{misc_line}\n')
        else:
            o.write(f'##fileformat=KCFv1.0 (Kmer CountTable Format)\n')
            o.write(f'##source=IBSpy\n')
            o.write(f'##TOTAL_KMER=Total number of Kmers counted within the window\n')
            o.write(f'##VA=Number of Variations\n')
            o.write(f'##OB=Observed k-mers\n')
            o.write(f'##DI=K-mer distance\n')
            o.write(f'##SC=Score calculated as (OB - VA)/TOTAL_KMER\n')
        o.write(f'##CMD: {" ".join(sys.argv)}\n')
        o.write(f'#CHROM\tSTART\tEND\tTOTAL_KMER\tINFO\tFORMAT\t{samples_line}\n')
        for window_key in windows:
            window = windows[window_key]
            o.write(f'{window}'
                    f'\t{windows[window_key].total_kmers}'
                    f'\t{windows[window_key].get_info()}'
                    f'\tIB:VA:OB:DI:SC'
                    f'\t{windows[window_key].get_window()}\n')


def write_bedgraph(windows, samples, output_prefix):
    """
    Write bedgraph files
    """
    for sample in samples:
        output_file = f'{output_prefix}.{sample}.bedgraph'
        logging.info(f'Writing {output_file}')
        fo = open(output_file, 'w')
        # fo.write(f'track type=bedGraph name="{sample}" description="{sample}" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n')
        prev_chrom = None
        for key in windows:
            chrom = key[0]
            start = key[1]
            end = key[2]
            value = windows[key].data[sample].score
            if prev_chrom != chrom:
                fo.write(f'{chrom}\t{start}\t{end}\t{value}\n')
            elif prev_end > start:
                fo.write(f'{chrom}\t{prev_end}\t{end}\t{value}\n')
            prev_end = end
            prev_chrom = chrom
        fo.close()


def write_matrix(windows, output, samples=None, allele_a_cutoff=0.8, allele_b_cutoff=0.3):
    """
    Write genotype matrix
    if the score is more than 0.8 make it 1,
    if the score is less than 0.3 make it 0,
    if the score is between 0.3 and 0.8 make it 0.5
    if the score is NA make it NA
    the output genotype matrix should have only the samples on the columns and window as rows
    """
    logging.info(f'Writing {output}.matrix.tr.tsv')
    with open(output + '.matrix.tr.tsv', 'w') as f_matrix, open(output + '.map.tsv', 'w') as f_map:
        f_matrix.write('taxa\t' + '\t'.join(samples) + '\n')
        f_map.write('name\tchromosome\tposition\n')
        for i, key in enumerate(windows):
            window = windows[key]
            scores = window.get_value('SC', samples)
            observed = window.get_value('OB', samples)
            genotypes = []
            for j, score in enumerate(scores):
                if score == 0.00:
                    if observed[j] == 0:
                        genotypes.append('N')
                    else:
                        genotypes.append('0')
                elif score >= allele_a_cutoff:
                    genotypes.append('2')
                elif score <= allele_b_cutoff:
                    genotypes.append('1')
                elif allele_b_cutoff < score < allele_a_cutoff:
                    genotypes.append('0')
                else:
                    genotypes.append('N')
            if 'N' not in genotypes and not monomorphic(genotypes):
                f_matrix.write(f'{i + 1}\t' + '\t'.join(genotypes) + '\n')
                f_map.write(f'{i + 1}\t{key[0]}\t{key[1]}\n')


def transpose_gt_matrix(input_file, output_file):
    """
    Transpose large tsv file
    """
    logging.info(f'Transposing {input_file}')
    # Get the total number of columns
    with open(input_file, 'r') as f:
        num_cols = len(f.readline().split('\t'))

    # Open the output file
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')

        for col_index in range(num_cols):
            with open(input_file, 'r') as infile:
                reader = csv.reader(infile, delimiter='\t')
                column = [row[col_index] for row in reader]

            # Write the column as a row in the output file
            writer.writerow(column)
