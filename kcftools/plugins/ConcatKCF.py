import logging
import os
import sys

from kcftools.utils.helper_functions import list_to_str
from kcftools._defaults import NEWLINE


class ConcatKCF:
    def __init__(self, args):
        self.kcfs = args.kcfs
        self.out_kcf = args.out_kcf

    def run(self):
        """
        Concatenate list of kcf files from different chromosomes to a single kcf file
        """
        fo = open(self.out_kcf, 'w')
        misc_lines = []
        first_file = True
        for in_kcf in self.kcfs:
            logging.info(f'Reading {in_kcf}')
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
                        fo.write(f'{list_to_str(misc_lines, NEWLINE)}\n')
                    else:
                        current_samples = line.strip().split('\t')[6:]
                        if samples != current_samples:
                            os.remove(self.out_kcf)
                            logging.info(f'Error: {in_kcf} contains different samples than previous files')
                            sys.exit(1)
                    continue
                fo.write(f'{line}\n')
            f.close()
        fo.close()