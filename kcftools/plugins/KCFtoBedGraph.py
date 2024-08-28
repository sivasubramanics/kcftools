#!/usr/bin/env python3

"""
Convert kcf file to bedgraph files
"""
import logging

from kcftools.utils.readers import read_kcf


class KCFtoBedGraph:
    def __init__(self, args):
        self.input_file = args.input
        self.output_prefix = args.output
        if args.sample:
            self.sample_name = args.sample
        else:
            self.sample_name = None

    def run(self):
        """
        Convert kcf file to bedgraph files
        """
        windows, samples, misc_lines = read_kcf(self.input_file)
        if self.sample_name is not None:
            if type(self.sample_name) == str:
                samples = [self.sample_name]
            elif type(self.sample_name) == list:
                samples = self.sample_name
        write_bedgraph(windows, samples, self.output_prefix)


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

