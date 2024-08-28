#!/usr/bin/env python3
# script to process IBSpy output
import os
import argparse
import logging
import subprocess
import sys
from kcftools._defaults import nthreads, mem

from kcftools.run import run


def main():
    """
    Main function
    """

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')

    # parse arguments for the "run_ibspy" command
    get_kmers = subparsers.add_parser('get_kmers', help='Run KMC for fasta or fastq files and generate a kmc database')
    get_kmers.add_argument('-i', '--input', help='fasta or fastq file(s)', type=str, nargs='+', required=False)
    get_kmers.add_argument('-I', '--input_list', help='file containing list of fasta or fastq files', type=str,
                           required=False)
    get_kmers.add_argument('-o', '--output', help='output directory', required=True)
    get_kmers.add_argument('-t', '--threads', help=f"number of threads [{nthreads}]", default=nthreads, type=int)
    get_kmers.add_argument('-m', '--mem', help=f"memory in GB [{mem}]", default=mem, type=int)
    get_kmers.add_argument('-f', '--format', help='format: fa - fasta single line. fq - fastq. fm - fasta '
                                                  'multiline', required=True, choices=['fa', 'fq', 'fm'])
    get_kmers.add_argument('-k', '--kmer', help='kmer size', required=True)

    run_ibspy = subparsers.add_parser('run_ibspy', help='Run IBSpy on the output of KMC and generate variant tables')
    run_ibspy.add_argument('-i', '--input', help='kmc database directory', required=True)
    run_ibspy.add_argument('-r', '--reference', help='reference fasta file', required=True)
    run_ibspy.add_argument('-n', '--name', help='reference name to be added with output', required=True)
    run_ibspy.add_argument('-k', '--kmer', help='kmer size', required=True)
    run_ibspy.add_argument('-w', '--window', help='window size', required=True)
    run_ibspy.add_argument('-t', '--threads', help=f"number of threads [{nthreads}]", default=nthreads, type=int)

    # Create the parser for the "tsv2vcf" command
    parser_tsv2vcf = subparsers.add_parser('tsv2kcf', help='Step 1: Convert IBSpy tsv to kcf')
    parser_tsv2vcf.add_argument('-i', '--input', help='Input IBSpy tsv file', required=True)
    parser_tsv2vcf.add_argument('-o', '--output', help='Output kcf file', required=True)
    parser_tsv2vcf.add_argument('-s', '--sample',
                                help='Sample name (if not given, will be taken from input file name)')

    # Create the parser for the "cohort" command
    parser_combine = subparsers.add_parser('cohort', help='Step 2.a: Combine kcf files')
    parser_combine.add_argument('-i', '--input', help='Input file containing list of VCF files', required=True)
    parser_combine.add_argument('-o', '--output', help='Output vcf file', required=True)

    # Create the parser for the "increase_window" command
    parser_increase_window = subparsers.add_parser('increase_window', help='Step 2.b: Increase window size')
    parser_increase_window.add_argument('-i', '--input', help='Input kcf file', required=True)
    parser_increase_window.add_argument('-o', '--output', help='Output kcf file', required=True)
    parser_increase_window.add_argument('-l', '--length',
                                        help='Sequence length file. Preferrably .fai', required=True)
    parser_increase_window.add_argument('-w', '--window', help='Window size', required=True)
    parser_increase_window.add_argument('-k', '--kmer',
                                        help='Kmer size used during the IBSpy pipeline', required=True)

    # Create the parser for the "find_IBS" command
    parser_find_IBS = subparsers.add_parser('find_IBS', help='Step 3: Find the IBS regions in kcf file')
    parser_find_IBS.add_argument('-i', '--input', help='Input kcf file', required=True)
    parser_find_IBS.add_argument('-o', '--output', help='Output kcf file', required=True)
    parser_find_IBS.add_argument('-v', '--variations', help='Minimum variations cut-off. [default: 6] ', default=6,
                                 type=float)
    parser_find_IBS.add_argument('-s', '--score', help='Minimum score cut-off. [default: 0.8]', default=0.8, type=float)
    parser_find_IBS.add_argument('-c', '--consecutive',
                                 help='Minimum consecutive windows with NA\'s [default: 5]', default=5, type=int)
    parser_find_IBS.add_argument('-r', '--reverse', help='set this if the IBS to be detected on contrast',
                                 action='store_true', default=False)

    # Create the parser for the "extract" command
    parser_extract = subparsers.add_parser('extract', help='Step 4: Extract windows from kcf file')
    parser_extract.add_argument('-i', '--input', help='Input kcf file', required=True)
    parser_extract.add_argument('-o', '--output', help='Output prefix', required=True)
    parser_extract.add_argument('-s', '--sample',
                                help='Sample name (if not given, will be taken from input file name)')
    parser_extract.add_argument('-H', '--heatmap', help='Output heatmap file', action='store_true', default=False)

    # Create the parser for the "bedgraph" command
    parser_bedgraph = subparsers.add_parser('kcf2bedgraph', help='Step 5: Convert kcf file to bedgraph files')
    parser_bedgraph.add_argument('-i', '--input', help='Input kcf file', required=True)
    parser_bedgraph.add_argument('-o', '--output', help='Output prefix', required=True)
    parser_bedgraph.add_argument('-s', '--sample',
                                 help='Sample name (if not given, will be taken from input file name)')

    # Create the parser for the "matrix" command
    parser_kcf2matrix = subparsers.add_parser('kcf2matrix', help='Step 6: Convert kcf file to genotype matrix')
    parser_kcf2matrix.add_argument('-i', '--input', help='Input kcf file', required=True)
    parser_kcf2matrix.add_argument('-o', '--output', help='Output prefix', required=True)
    parser_kcf2matrix.add_argument('-s', '--sample',
                                   default=None,
                                   help='Sample name (if not given, will be taken from input file name)')
    parser_kcf2matrix.add_argument('-a', '--allele_a_cutoff', help='Allele A cutoff [default: 0.8]', default=0.8,
                                   type=float)
    parser_kcf2matrix.add_argument('-b', '--allele_b_cutoff', help='Allele B cutoff [default: 0.2]', default=0.2,
                                   type=float)

    # Create the parser for the "split_kcf" command
    parser_split_kcf = subparsers.add_parser('split_kcf', help='Misc: Split kcf file by chromosome')
    parser_split_kcf.add_argument('-i', '--input', help='Input kcf file', required=True)
    parser_split_kcf.add_argument('-o', '--output', help='Output prefix', required=True)
    parser_split_kcf.add_argument('-s', '--sample',
                                  default=None,
                                  help='Sample name (if not given, will be taken from input file name)')
    parser_split_kcf.add_argument('-c', '--chrs', help='Chromosomes to be extracted', nargs='+')

    # Create the parser for the "concat" command
    parser_concat = subparsers.add_parser('concat', help='Misc: Concatenate kcf files (samples should be identical)')
    parser_concat.add_argument('-i', '--input', help='Input kcf files', nargs='+', required=True)
    parser_concat.add_argument('-o', '--output', help='Output kcf file', required=True)

    # Parser for extract_samples command
    parser_extract_samples = subparsers.add_parser('extract_samples', help='Misc: Extract samples from kcf file')
    parser_extract_samples.add_argument('-i', '--input', help='Input kcf file (multi sampled)', required=True)
    parser_extract_samples.add_argument('-o', '--output', help='Output kcf file', required=True)
    parser_extract_samples.add_argument('-s', '--samples',
                                        help='Samples to be extracted Names or list',
                                        nargs='+', default=None)
    parser_extract_samples.add_argument('-S', '--samples_list', help='File containing list of samples to be extracted',
                                        default=None)

    parser_get_scores = subparsers.add_parser('get_attr', help='Misc: Get values of the attributes from kcf file')
    parser_get_scores.add_argument('-i', '--input', help='Input kcf file', required=True)
    parser_get_scores.add_argument('-o', '--output', help='Output file', required=True)
    parser_get_scores.add_argument('-a', '--attribute', help='Attribute to be extracted', required=True,
                                   choices=['VA', 'OB', 'DI', 'SC'])
    parser_get_scores.add_argument('-s', '--sample',
                                   help='Sample name (if not given, will be taken from input file name)', nargs='+',
                                   default=None)
    parser_get_scores.add_argument('-S', '--samples_list', help='File containing list of samples to be extracted',
                                   default=None)

    parser_extract_score_matrix = subparsers.add_parser('extract_score_matrix',
                                                        help='Misc: Extract score matrix from ibs.tsv files (output of extract step)')
    parser_extract_score_matrix.add_argument('-i', '--input', help='TSV file with label and ibs.tsv file',
                                             required=True)
    parser_extract_score_matrix.add_argument('-o', '--output', help='Output file', required=True)
    parser_extract_score_matrix.add_argument('-l', '--length', help='genome length file (.fai file)', required=True)
    parser_extract_score_matrix.add_argument('-m', '--min_len', help='minimum length to filter', default=1000000,
                                             type=int)
    parser_extract_score_matrix.add_argument('-s', '--score', help='minimum score to filter. [default: 0.8]',
                                             default=0.8, type=float)
    parser_extract_score_matrix.add_argument('-r', '--reverse', help='set this if the IBS to be detected on contrast',
                                             action='store_true', default=False)
    parser_extract_score_matrix.add_argument('-p', '--ibs_prop',
                                             help='minimum IBS blocks proportion to filter. [default: 0.5]',
                                             default=0.5, type=float)

    parser_count_genes = subparsers.add_parser('count_genes', help='Misc: Count genes in the IBS blocks')
    parser_count_genes.add_argument('-i', '--input', help='Input TSV file from extract step', required=True)
    parser_count_genes.add_argument('-o', '--output', help='Output file', required=True)
    parser_count_genes.add_argument('-g', '--gtf', help='GTF file', required=True)
    parser_count_genes.add_argument('-f', '--feature', help='Feature to be counted', default='gene')
    parser_count_genes.add_argument('-s', '--score', help='minimum score to filter. [default: 0.8]', default=0.8,
                                    type=float)
    parser_count_genes.add_argument('-r', '--reverse', help='set this if the IBS to be detected on contrast',
                                    action='store_true', default=False)
    parser_count_genes.add_argument('-p', '--ibs_prop', help='minimum IBS blocks proportion to filter. [default: 0.5]',
                                    default=0.5, type=float)
    parser_count_genes.add_argument('-m', '--min_length', help='minimum length to filter', default=10000, type=int)

    parse_recalculate = subparsers.add_parser('recalculate', help='Misc: Recalculate the scores')
    parse_recalculate.add_argument('-i', '--input', help='Input kcf file', required=True)
    parse_recalculate.add_argument('-o', '--output', help='Output kcf file', required=True)

    args = parser.parse_args(args=(sys.argv[1:] or ['--help']))

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    run(args)
