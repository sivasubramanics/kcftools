#!/usr/bin/env python3

"""
Find IBS regions
"""

import logging

from kcftools.utils.readers import read_kcf
from kcftools.utils.writers import write_kcf


class FindIBS:
    def __init__(self, args):
        self.input_file = args.input
        self.output_file = args.output
        self.min_variations = args.variations
        self.min_score = args.score
        self.min_consecutive = args.consecutive
        self.reverse = args.reverse

    def run(self):
        """
        Flag the IBS regions in kcf file
        """
        windows, samples, misc_lines = read_kcf(self.input_file)
        logging.info(f'Finding IBS regions in {self.input_file}')
        for sample in samples:
            block_num = 0
            na_num = 0
            last_chrom = None
            first_block = True
            for key in windows:
                if not self.reverse:
                    if key[0] != last_chrom and last_chrom is not None:
                        first_block = True
                        block_num += 1
                        na_num = 0
                    if windows[key].data[sample].score >= self.min_score:
                        if first_block:
                            block_num += 1
                        if na_num > self.min_consecutive:
                            block_num += 1
                        na_num = 0
                        windows[key].data[sample].ibs = str(block_num)
                    elif windows[key].data[sample].va <= self.min_variations and windows[key].data[sample].score >= 0.5:
                        if first_block:
                            block_num += 1
                        if na_num > self.min_consecutive:
                            block_num += 1
                        na_num = 0
                        windows[key].data[sample].ibs = str(block_num)
                    else:
                        windows[key].data[sample].ibs = 'N'
                        na_num += 1
                    last_chrom = key[0]
                    first_block = False
                else:
                    if key[0] != last_chrom and last_chrom is not None:
                        first_block = True
                        block_num += 1
                        na_num = 0
                    if windows[key].data[sample].score <= self.min_score and windows[key].data[sample].score >= 0.5:
                        if first_block:
                            block_num += 1
                        if na_num > self.min_consecutive:
                            block_num += 1
                        na_num = 0
                        windows[key].data[sample].ibs = str(block_num)
                    elif windows[key].data[sample].va >= self.min_variations and windows[key].data[sample].score <= 0.2:
                        if first_block:
                            block_num += 1
                        if na_num > self.min_consecutive:
                            block_num += 1
                        na_num = 0
                        windows[key].data[sample].ibs = str(block_num)
                    else:
                        windows[key].data[sample].ibs = 'N'
                        na_num += 1
                    last_chrom = key[0]
                    first_block = False
        write_kcf(self.output_file, windows, samples, misc_lines)
