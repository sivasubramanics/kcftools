#!/usr/bin/env python3
"""
Extract samples from kcf file
"""

import sys
import logging

from kcftools.utils.helper_functions import list_to_str
from kcftools._defaults import TAB
from kcftools.utils.parsers import update_info


class ExtractSamples:
    def __init__(self, args):
        self.input = args.input
        self.output = args.output
        self.samples = args.samples

    def run(self):
        """
        Extract samples from kcf file
        """
        # TODO: fix this function for modifying the info fields as well
        sample_indices = []
        if type(self.samples) == str:
            samples = [self.samples]
        elif type(self.samples) == list:
            samples = self.samples
        else:
            logging.error('Error: samples should be a string or a list of strings')
            sys.exit(1)
        fi = open(self.input, 'r')
        fo = open(self.output, 'w')
        for line in fi:
            if line.startswith('##'):
                fo.write(line)
                continue
            line = line.strip().split('\t')
            if line[0].startswith('#CHROM'):
                fo.write(f"##CMD: {' '.join(sys.argv)}\n")
                for sample in samples:
                    if sample in line:
                        sample_indices.append(line.index(sample))
                fo.write(f'{list_to_str(line[:6], TAB)}')
                for i in sample_indices:
                    fo.write(f'\t{line[i]}')
                fo.write('\n')
                continue
            line[4] = update_info(line, sample_indices)
            fo.write(f'{list_to_str(line[:6])}\t{list_to_str([line[i] for i in sample_indices])}\n')
        fi.close()
        fo.close()