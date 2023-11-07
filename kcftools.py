#!/usr/bin/env python3

import os
import sys
import argparse
import time
from collections import defaultdict, Counter
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
import csv


NEWLINE = '\n'
TAB = '\t'
COMMA = ','
COLON = ':'
SEMICOLON = ';'

def print_log(message):
    """
    Print log message
    """
    print(f'[{time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}] {message}')
    sys.stdout.flush()


class Data:
    def __init__(self, ibs, ob, va, di, tot):
        self.ob = ob
        self.va = va
        self.di = di
        self.tot = tot
        self.ibs = ibs

    @property
    def score(self):
        if self.tot == 0:
            return 0.00
        if self.ob == 0:
            return 0.00
        _score = round((self.ob - self.va) / self.tot, 2)
        if _score < 0.001:
            return 0.00
        return _score

    def __str__(self):
        return f'{self.ibs}:{self.va}:{self.ob}:{self.di}:{self.score}'


class Window:
    def __init__(self, chr_name, start, end):
        self.chr_name = chr_name
        self.start = start
        self.end = end
        self.total_kmers = 0
        self.data = defaultdict()
        self.blocks = defaultdict()

    def set_total_kmers(self, total_kmers):
        self.total_kmers = total_kmers

    def add_total_kmers(self, total_kmers):
        self.total_kmers += total_kmers

    def set_value(self, sample_names, attribute, value):
        if type(sample_names) == str:
            sample_names = [sample_names]
        for sample_name in sample_names:
            if sample_name not in self.data:
                self.data[sample_name] = Data(0, 0, 0, self.total_kmers)
            if attribute == 'VA':
                self.data[sample_name].va = value
            elif attribute == 'OB':
                self.data[sample_name].ob = value
            elif attribute == 'DI':
                self.data[sample_name].di = value

    def set_data(self, sample_name, data):
        if sample_name not in self.data:
            self.data[sample_name] = Data('N', 0, 0, 0, self.total_kmers)
        self.data[sample_name] = data

    def add_data(self, sample_name, data):
        if sample_name not in self.data:
            self.data[sample_name] = Data('N', 0, 0, 0, 0)
        self.data[sample_name].tot += data.tot
        self.data[sample_name].ob += data.ob
        self.data[sample_name].va += data.va
        self.data[sample_name].di += data.di

    def get_value(self, attribute, sample_names=None):
        if sample_names is None:
            sample_names = self.data.keys()
        if type(sample_names) == str:
            sample_names = [sample_names]
        values = []
        for sample_name in sample_names:
            if attribute == 'VA':
                values.append(self.data[sample_name].va)
            elif attribute == 'OB':
                values.append(self.data[sample_name].ob)
            elif attribute == 'DI':
                values.append(self.data[sample_name].di)
            elif attribute == 'SC':
                values.append(self.data[sample_name].score)
            elif attribute == 'IB':
                values.append(self.data[sample_name].ibs)
            else:
                sys.exit(f'Error: Attribute {attribute} is not valid')
        return values

    def get_window(self, sample_names=None):
        if sample_names is None:
            sample_names = self.data.keys()
        vcf_list = []
        if type(sample_names) == str:
            sample_names = [sample_names]
        for sample_name in sample_names:
            vcf_list.append(f"{self.data[sample_name]}")
        return '\t'.join(vcf_list)

    def get_info(self, sample_names=None):
        if sample_names is None:
            sample_names = self.data.keys()
        if type(sample_names) == str:
            sample_names = [sample_names]
        scores = self.get_value('SC', sample_names)
        va = self.get_value('VA', sample_names)
        ob = self.get_value('OB', sample_names)
        di = self.get_value('DI', sample_names)
        return (f"MIN_SCORE={min(scores)};"
                f"MAX_SCORE={max(scores)};"
                f"MEAN_SCORE={round(sum(scores) / len(scores), 2)};"
                f"MIN_VA={min(va)};"
                f"MAX_VA={max(va)};"
                f"MEAN_VA={round(sum(va) / len(va), 2)};"
                f"MIN_OB={min(ob)};"
                f"MAX_OB={max(ob)};"
                f"MEAN_OB={round(sum(ob) / len(ob), 2)};")

    def __str__(self):
        return f'{self.chr_name}\t{self.start}\t{self.end}'


def duplicates(it):
    """
    Check if there are duplicates in a list and return the duplicates
    """
    return [item for item, count in Counter(it).items() if count > 1]


def read_tsv(input_file, sample_name=None):
    """
    Read IBSpy tsv file and return a Window object
    """
    if sample_name is None:
        sample_name = os.path.basename(input_file).split('.')[0]
    line_count = 0
    print_log(f'Reading {input_file}')
    with open(input_file, 'r') as f:
        for line in f:
            line_count += 1
            if line_count == 1:
                if line.startswith('seqname'):
                    continue
                else:
                    sys.exit(f'Error: {input_file} does not appear to be a valid IBSpy tsv file')
            line = line.strip().split('\t')
            window = Window(line[0], int(line[1]), int(line[2]))
            window.set_total_kmers(int(line[3]))
            window.set_data(sample_name, Data('N', int(line[4]), int(line[5]), int(line[6]), int(line[3])))
            yield window


def tsv2kcf(input_file, output_file, sample_name=None):
    """
    Convert IBSpy tsv file to kcf file
    """
    if sample_name is None:
        sample_name = os.path.basename(input_file).split('.')[0]
    windows = defaultdict()
    with open(output_file, 'w') as o:
        for window in read_tsv(input_file, sample_name):
            windows[(window.chr_name, window.start, window.end)] = window
    write_kcf(output_file, windows, sample_name)


def get_chromosomes(windows):
    """
    Return a list of chromosomes in the windows
    """
    chromosomes = []
    for window in windows:
        if window[0] not in chromosomes:
            chromosomes.append(window[0])
    return chromosomes


def read_kcf(input_file, windows=None, samples_list=None):
    """
    Read kcf file and return a Window object
    """
    if windows is None:
        windows = defaultdict(Window)
    if samples_list is None:
        samples_list = []
    line_count = 0
    misc_lines = []
    print_log(f'Reading {input_file}')
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('##'):
                misc_lines.append(line.strip())
                continue
            line_count += 1
            if line_count == 1:
                if line.startswith('#CHROM'):
                    samples = line.strip().split('\t')[6:]
                    if duplicates(samples):
                        sys.exit(f'Error: {input_file} contains duplicate samples')
                    # samples_list.extend(samples)
                    samples_list += samples
                    if duplicates(samples_list):
                        sys.exit(f'Error: {input_file} contains samples that are already present in other input files')
                    continue
                else:
                    sys.exit(f'Error: {input_file} does not appear to be a valid kcf file')
            line = line.strip().split('\t')
            key = (line[0], int(line[1]), int(line[2]))
            window = Window(line[0], int(line[1]), int(line[2]))
            if key not in windows:
                windows[key] = window
            windows[key].set_total_kmers(int(line[3]))
            for i, sample in enumerate(samples):
                ib, va, ob, di, sc = line[6 + i].split(':')
                if ib != 'N':
                    ib = int(ib)
                windows[key].set_data(sample, Data(ib, int(ob), int(va), int(di), int(line[3])))
    return windows, samples_list, misc_lines


def cohort(input_files, output_file):
    """
    Combine multiple kcf files
    """
    windows = defaultdict()
    samples = []
    with open(input_files, 'r') as f:
        for input_file in f:
            input_file = input_file.strip()
            windows, samples, misc_lines = read_kcf(input_file, windows, samples)
    write_kcf(output_file, windows, samples, misc_lines)


def increase_window(input_file, output_file, window_size, kmer_size, length_file):
    """
    Increase window size
    """
    windows, samples, misc_lines = read_kcf(input_file)
    seq_lengths = get_seq_lengths(length_file)
    current_window_size = get_current_window_size(windows, seq_lengths)
    print_log(f'Current window size: {current_window_size}')
    windows = convert_windows(windows, seq_lengths, current_window_size, int(window_size), int(kmer_size))
    print_log(f'New window size: {window_size}')
    write_kcf(output_file, windows, samples, misc_lines)


def write_kcf(out_kcf, windows, samples, misc_lines=None):
    """
    Write kcf file
    """
    if type(samples) == str:
        samples = [samples]
    print_log(f'Writing {out_kcf}')
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


def get_current_window_size(windows, chr_lengths):
    """
    Get current window size
    """
    window_size = 0
    for (seqname, start, end) in windows:
        if window_size == 0:
            window_size = end - start
        if end == chr_lengths[seqname]:
            continue
        if window_size != end - start:
            sys.exit(f'Error: Window size is not consistent in input file.')
    return int(window_size)


def convert_windows(windows, chr_lengths, o_win_size, n_win_size, kmer_size):
    """
    Convert window size
    """
    # get all the chromosomes from the windows dictionary
    chromosomes = get_chromosomes(windows)
    new_windows = defaultdict(dict)
    for seqname in chr_lengths:
        # check if the chromosome is in the fai file
        if seqname not in chromosomes:
            continue
        n_win_start = 0
        win_starts, win_ends = get_window(0, chr_lengths[seqname], n_win_size, kmer_size)
        for win_start, win_end in zip(win_starts, win_ends):
            new_window = Window(seqname, win_start, win_end)
            new_key = (str(seqname), str(win_start), str(win_end))
            o_win_starts, o_win_ends = get_window(n_win_start, win_end, o_win_size, kmer_size)
            for o_win_start, o_win_end in zip(o_win_starts, o_win_ends):
                key = (str(seqname), o_win_start, o_win_end)
                if key not in windows:
                    n_win_start = o_win_start
                    continue
                old_window = windows[key]
                samples = old_window.data.keys()
                new_window.add_total_kmers(old_window.total_kmers)
                for sample in samples:
                    new_window.add_data(sample, old_window.data[sample])
            new_windows[new_key] = new_window
    return new_windows


def get_window(start, end, window_size, kmer_size):
    """
    return a list of windows
    """
    win_start = []
    win_end = []
    t_win_start = start
    t_win_end = start + window_size
    while t_win_end < end:
        win_start.append(t_win_start)
        win_end.append(t_win_end)
        t_win_start = t_win_end - kmer_size
        t_win_end = t_win_start + window_size
    if t_win_end >= end:
        win_start.append(t_win_start)
        win_end.append(end)
    return win_start, win_end


def get_seq_lengths(length_file):
    """
    Read sequence length file and return a dictionary
    """
    chr_lengths = defaultdict()
    with open(length_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            chr_lengths[line[0]] = int(line[1])
    return chr_lengths


def find_IBS(input_file, output_file, min_variations, min_score, min_consecutive):
    """
    Flag the IBS regions in kcf file
    """
    windows, samples, misc_lines = read_kcf(input_file)
    print_log(f'Finding IBS regions in {input_file}')
    for sample in samples:
        block_num = 0
        na_num = 0
        last_chrom = None
        first_block = True
        for key in windows:
            if key[0] != last_chrom and last_chrom is not None:
                first_block = True
                block_num += 1
                na_num = 0
            if windows[key].data[sample].score >= min_score:
                if first_block:
                    block_num += 1
                if na_num > min_consecutive:
                    block_num += 1
                na_num = 0
                windows[key].data[sample].ibs = str(block_num)
            elif windows[key].data[sample].va <= min_variations and windows[key].data[sample].score >= 0.5:
                if first_block:
                    block_num += 1
                if na_num > min_consecutive:
                    block_num += 1
                na_num = 0
                windows[key].data[sample].ibs = str(block_num)
            else:
                windows[key].data[sample].ibs = 'N'
                na_num += 1
            last_chrom = key[0]
            first_block = False
    write_kcf(output_file, windows, samples, misc_lines)


def extract(input_file, output_prefix, sample_name=None):
    """
    Extract windows from kcf file
    """
    windows, samples, misc_lines = read_kcf(input_file)
    print_log(f'Extracting windows from {input_file}')
    if sample_name is not None:
        if type(sample_name) == str:
            samples = [sample_name]
        elif type(sample_name) == list:
            samples = sample_name

    for sample in samples:
        output_file = f'{output_prefix}.{sample}.tsv'
        fo = open(output_file, 'w')
        fo.write(f'seqname\tstart\tend\tlength\ttotal_blocks\tibs_blocks\n')
        out_bed = f'{output_prefix}.{sample}.bed'
        fo_bed = open(out_bed, 'w')
        last_block_num = 0
        block_count = 0
        block_start = None
        na_num = 0
        all_num = 0
        for (seqname, start, end) in windows:
            block_num = windows[(seqname, start, end)].data[sample].ibs
            all_num += 1
            if block_num == 'N':
                na_num += 1
                continue
            if block_num != last_block_num:
                if block_start is not None:
                    fo.write(
                        f'{seqname}\t{block_start}\t{block_end}\t{block_end - block_start}\t{all_num - na_num - 1}\t{block_count}\n')
                    fo_bed.write(f'{seqname}\t{block_start}\t{block_end}\t0\t+\n')
                block_count = 0
                block_start = start
                all_num = 1
            na_num = 0
            last_block_num = block_num
            block_end = end
            block_count += 1
        fo.write(
            f'{seqname}\t{block_start}\t{block_end}\t{block_end - block_start}\t{all_num - na_num}\t{block_count}\n')
        fo_bed.write(f'{seqname}\t{block_start}\t{block_end}\t0\t+\n')
        fo.close()

    prev_chrom = None
    num_window = 0
    for key in windows:
        chrom = key[0]
        num_window += 1
        if chrom != prev_chrom:
            if prev_chrom is not None:
                fo.close()
            num_window = 1
            fo = open(f'{output_prefix}.{chrom}.heatmap.tsv', 'w')
            samples_line = '\t'.join(samples)
            fo.write(f'window\t{samples_line}\n')
        ibs = "\t".join(map(ibs_to_binary, windows[key].get_value('IB', samples)))
        fo.write(f'{num_window}\t{ibs}\n')
        prev_chrom = chrom
    fo.close()


def kcf2bedgraph(input_file, output_prefix, sample_name=None):
    """
    Convert kcf file to bedgraph files
    """
    windows, samples, misc_lines = read_kcf(input_file)
    if sample_name is not None:
        if type(sample_name) == str:
            samples = [sample_name]
        elif type(sample_name) == list:
            samples = sample_name
    write_bedgraph(windows, samples, output_prefix)


def write_bedgraph(windows, samples, output_prefix):
    """
    Write bedgraph files
    """
    for sample in samples:
        output_file = f'{output_prefix}.{sample}.bedgraph'
        print_log(f'Writing {output_file}')
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


def ibs_to_binary(in_str):
    """
    Convert IBS to binary
    """
    if in_str == 'N':
        return '0'
    else:
        return '1'


def write_matrix(windows, output, samples=None):
    """
    Write genotype matrix
    if the score is more than 0.8 make it 1,
    if the score is less than 0.3 make it 0,
    if the score is between 0.3 and 0.8 make it 0.5
    if the score is NA make it NA
    the output genotype matrix should have only the samples on the columns and window as rows
    """
    print_log(f'Writing {output}.matrix.tr.tsv')
    with open(output + '.matrix.tr.tsv', 'w') as f_matrix, open(output + '.map.tsv', 'w') as f_map:
        f_matrix.write('taxa\t' + '\t'.join(samples) + '\n')
        f_map.write('name\tchromosome\tposition\n')
        for i, key in enumerate(windows):
            window = windows[key]
            scores = window.get_value('SC', samples)
            observed = window.get_value('OB', samples)
            genotypes = []
            for i, score in enumerate(scores):
                if score == 0.00:
                    if observed[i] == 0:
                        genotypes.append('N')
                    else:
                        genotypes.append('0')
                elif score >= 0.80:
                    genotypes.append('2')
                elif score <= 0.30:
                    genotypes.append('0')
                elif 0.30 < score < 0.80:
                    genotypes.append('1')
                else:
                    genotypes.append('N')
            f_matrix.write(f'{i + 1}\t' + '\t'.join(genotypes) + '\n')
            f_map.write(f'{i+1}\t{key[0]}\t{key[1]}\n')


def transpose_gt_matrix(input_file, output_file):
    """
    Transpose large tsv file
    """
    print_log(f'Transposing {input_file}')
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

def kcf2matrix(ikcf, outprefix, sample):
    """
    Convert kcf file to genotype matrix and genotype map file.
    the map file will have the window name, chromosome and the position in TSV format
    """
    windows, samples, misc_lines = read_kcf(ikcf)
    if sample is not None:
        if type(sample) == str:
            samples = [sample]
        elif type(sample) == list:
            samples = sample
    write_matrix(windows, outprefix, samples)
    transpose_gt_matrix(f"{outprefix}.matrix.tr.tsv", f"{outprefix}.matrix.tsv")


def list_to_str(in_list, sep='\t'):
    """
    Convert list to string
    """
    return sep.join(map(str, in_list))


def split_kcf(in_kcf, out_prefix, samples=None, chrs=None):
    """
    Split kcf file by chromosome for given set of samples if list of samples provided
    """
    prev_chrom = None
    data_indices = [0, 1, 2, 3, 4, 5]
    misc_lines = []
    with open(in_kcf, 'r') as f:
        for line in f:
            if line.startswith('##'):
                misc_lines.append(line.strip())
                continue
            if line.startswith('#CHROM'):
                misc_lines.append(f"##CMD: {' '.join(sys.argv)}")
                in_samples = line.strip().split('\t')[6:]
                if samples is None:
                    samples = in_samples
                if type(samples) == str:
                    samples = [samples]
                for sample in samples:
                    data_indices.append(in_samples.index(sample)+6)
                    if sample not in in_samples:
                        sys.exit(f'Error: {sample} not found in {in_kcf}')
                misc_lines.append(f'#CHROM\tSTART\tEND\tTOTAL_KMER\tINFO\tFORMAT\t{list_to_str(samples)}')
                continue
            else:
                line = line.strip().split('\t')
                out_line = [line[i] for i in data_indices]
                if chrs is not None:
                    if line[0] not in chrs:
                        continue
                if prev_chrom is None or line[0] != prev_chrom:
                    print_log(f'Writing {out_prefix}.{line[0]}.kcf')
                    fo = open(f'{out_prefix}.{line[0]}.kcf', 'w')
                    fo.write(f'{list_to_str(misc_lines, NEWLINE)}\n')
                fo.write(f'{list_to_str(out_line)}\n')
                prev_chrom = line[0]


def concat(in_kcfs, out_kcf):
    """
    Concatenate list of kcf files from different chromosomes to a single kcf file
    """
    fo = open(out_kcf, 'w')
    misc_lines = []
    first_file = True
    for in_kcf in in_kcfs:
        print_log(f'Reading {in_kcf}')
        f = open(in_kcf, 'r')
        for line in f:
            line = line.strip()
            if line.startswith('##'):
                misc_lines.append(line)
                continue
            if line.startswith('#CHROM'):
                if first_file:
                    samples = line.strip().split('\t')[6:]
                    misc_lines.append(f'#CHROM\tSTART\tEND\tTOTAL_KMER\tINFO\tFORMAT\t{list_to_str(samples)}')
                    first_file = False
                else:
                    current_samples = line.strip().split('\t')[6:]
                    if samples != current_samples:
                        os.remove(out_kcf)
                        print_log(f'Error: {in_kcf} contains different samples than previous files')
                        sys.exit(1)
                fo.write(f'{list_to_str(misc_lines, NEWLINE)}\n')
                continue
            fo.write(f'{line}\n')
        f.close()
    fo.close()




def main():
    # start the clock
    start_time = time.time()

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')

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
    parser_find_IBS.add_argument('-v', '--variations', help='Minimum variations cut-off', default=6, type=float)
    parser_find_IBS.add_argument('-s', '--score', help='Minimum score cut-off', default=0.8, type=float)
    parser_find_IBS.add_argument('-c', '--consecutive',
                                 help='Minimum consecutive windows with NA\'s', default=5, type=int)

    # Create the parser for the "extract" command
    parser_extract = subparsers.add_parser('extract', help='Step 4: Extract windows from kcf file')
    parser_extract.add_argument('-i', '--input', help='Input kcf file', required=True)
    parser_extract.add_argument('-o', '--output', help='Output prefix', required=True)
    parser_extract.add_argument('-s', '--sample',
                                help='Sample name (if not given, will be taken from input file name)')

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
                                   help='Sample name (if not given, will be taken from input file name)')

    # Create the parser for the "split_kcf" command
    parser_split_kcf = subparsers.add_parser('split_kcf', help='Misc: Split kcf file by chromosome')
    parser_split_kcf.add_argument('-i', '--input', help='Input kcf file', required=True)
    parser_split_kcf.add_argument('-o', '--output', help='Output prefix', required=True)
    parser_split_kcf.add_argument('-s', '--sample',
                                  help='Sample name (if not given, will be taken from input file name)')
    parser_split_kcf.add_argument('-c', '--chrs', help='Chromosomes to be extracted', nargs='+')

    # Create the parser for the "concat" command
    parser_concat = subparsers.add_parser('concat', help='Misc: Concatenate kcf files (samples should be identical)')
    parser_concat.add_argument('-i', '--input', help='Input kcf files', nargs='+', required=True)
    parser_concat.add_argument('-o', '--output', help='Output kcf file', required=True)

    args = parser.parse_args(args=(sys.argv[1:] or ['--help']))

    if args.command == 'tsv2kcf':
        tsv2kcf(args.input, args.output, args.sample)
    elif args.command == 'cohort':
        cohort(args.input, args.output)
    elif args.command == 'increase_window':
        increase_window(args.input, args.output, args.window, args.kmer, args.length)
    elif args.command == 'find_IBS':
        find_IBS(args.input, args.output, args.variations, args.score, args.consecutive)
    elif args.command == 'extract':
        extract(args.input, args.output, args.sample)
    elif args.command == 'kcf2bedgraph':
        kcf2bedgraph(args.input, args.output, args.sample)
    elif args.command == 'kcf2matrix':
        kcf2matrix(args.input, args.output, args.sample)
    elif args.command == 'split_kcf':
        split_kcf(args.input, args.output, args.sample, args.chrs)
    elif args.command == 'concat':
        concat(args.input, args.output)
    else:
        parser.print_help()
        sys.exit(1)

    # stop clock
    end_time = time.time()
    print(f'Total time taken: {round(end_time - start_time, 2)} seconds')


if __name__ == '__main__':
    main()

