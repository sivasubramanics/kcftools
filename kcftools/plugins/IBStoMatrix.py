#!/usr/bin/env python3

"""
convert IBS output to matrix for heatmap purpose
"""

import logging
import os
from collections import defaultdict


class IBStoMatrix:
    def __init__(self, args):
        self.input = args.input
        self.output = args.output
        self.length = args.length
        self.min_length = args.min_length
        self.score_cutoff = args.score_cutoff
        self.ibs_prop = args.ibs_prop
        self.reverse = args.reverse

    def run(self):
        """
        Extract score matrix from IBS TSv files
        """
        len_dict = defaultdict()
        with open(self.length, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                len_dict[line[0]] = int(line[1])

        data_dict = defaultdict()
        with open(self.input, 'r') as f:
            for line in f:
                line = line.strip().split('\t')
                # check if file exits or accessible
                if not os.path.isfile(line[1]):
                    logging.info(f'Error: {line[1]} not found or not accessible. skipping')
                    continue
                data_dict[line[0]] = line[1]

        samples = list(data_dict.keys())
        logging.info(f'Total samples: {len(samples)}')

        out_matrix = defaultdict()
        for seqname in len_dict:
            for i in range(0, len_dict[seqname], 500000):
                win = int(i / 500000) * 0.5
                out_matrix[(seqname, win)] = [0] * len(samples)

        for sample in samples:
            logging.info(f'Processing {sample}')
            with open(data_dict[sample], 'r') as f:
                for line in f:
                    line = line.strip().split('\t')
                    if line[0] == 'id':
                        continue
                    if float(line[4]) < self.min_length:
                        continue
                    window_start = int(float(line[2]) / 500000)
                    window_end = int(float(line[3]) / 500000)
                    if float(line[6]) / float(line[5]) < self.ibs_prop:
                        continue
                    if self.reverse:
                        if float(line[7]) <= self.score_cutoff:
                            for i in range(window_start, window_end + 1):
                                out_matrix[(line[1], i * 0.5)][samples.index(sample)] = line[7]
                    else:
                        if float(line[7]) >= self.score_cutoff:
                            for i in range(window_start, window_end + 1):
                                out_matrix[(line[1], i * 0.5)][samples.index(sample)] = line[7]

        with open(self.output, 'w') as f:
            f.write('seqname\tstart\t' + '\t'.join(samples) + '\n')
            for key in out_matrix:
                f.write(f'{key[0]}\t{key[1]}\t' + '\t'.join(map(str, out_matrix[key])) + '\n')