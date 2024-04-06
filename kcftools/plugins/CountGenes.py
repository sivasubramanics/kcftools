#!/usr/bin/env python3

"""
Count the number of genes in each block
"""

import logging
from collections import defaultdict


class CountGenes:
    def __init__(self, args):
        self.input = args.input
        self.output = args.output
        self.gtf = args.gtf
        self.feature = args.feature
        self.min_length = args.min_length
        self.ibs_prop = args.ibs_prop
        self.reverse = args.reverse
        self.score_cutoff = args.score_cutoff

    def run(self):
        """
        Count the number genes per block
        """
        logging.info(f'Counting genes in {input}')
        genes = defaultdict()
        with open(self.gtf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip().split('\t')
                if line[0] not in genes:
                    genes[line[0]] = defaultdict()
                if line[2] == self.feature:
                    if line[3] not in genes[line[0]]:
                        genes[line[0]][line[3]] = []
                    genes[line[0]][line[3]].append(float(line[4]) - float(line[3]))

        with open(self.input, 'r') as fi, open(self.output, 'w') as fo:
            for line in fi:
                line = line.strip().split('\t')
                if line[2] == 'seqname':
                    line = '\t'.join(line)
                    fo.write(f"{line}\tn_genes\n")
                    continue
                n_genes = 0
                if float(line[4]) - float(line[3]) < self.min_length:
                    continue
                if float(line[7]) / float(line[6]) < self.ibs_prop:
                    continue
                if self.reverse:
                    if float(line[8]) > self.score_cutoff:
                        continue
                else:
                    if float(line[8]) < self.score_cutoff:
                        continue

                if line[2] in genes:
                    for g_start in genes[line[2]]:
                        for len in genes[line[2]][g_start]:
                            if float(g_start) >= float(line[3]) and float(g_start) + len <= float(line[4]):
                                n_genes += 1
                line = '\t'.join(line)
                fo.write(f"{line}\t{n_genes}\n")
