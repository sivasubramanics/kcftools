#!/usr/bin/env python3

"""
Readers for different file formats
"""

from kcftools.data.Window import Window, Data
from kcftools.utils.helper_functions import duplicates
import logging
import os
import sys
from collections import defaultdict


def read_tsv(input_file, sample_name=None):
    """
    Read IBSpy tsv file and return a Window object
    """
    if sample_name is None:
        sample_name = os.path.basename(input_file).split('.')[0]
    line_count = 0
    logging.info(f'Reading {input_file}')
    # check if the tsv file is empty
    if os.stat(input_file).st_size == 0:
        logging.error(f'Error: {input_file} is empty')
        sys.exit(1)
    with open(input_file, 'r') as f:
        for line in f:
            line_count += 1
            if line_count == 1:
                if line.startswith('seqname'):
                    continue
                else:
                    logging.error(f'Error: {input_file} does not appear to be a valid IBSpy tsv file')
                    sys.exit(1)
            line = line.strip().split('\t')
            window = Window(line[0], int(line[1]), int(line[2]))
            window.set_total_kmers(int(line[3]))
            window.set_data(sample_name, Data('N', int(line[4]), int(line[5]), int(line[6]), int(line[3])))
            yield window


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
    logging.info(f'Reading {input_file}')
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
