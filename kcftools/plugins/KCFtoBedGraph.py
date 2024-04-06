#!/usr/bin/env python3

"""
Convert kcf file to bedgraph files
"""

from kcftools.utils.readers import read_kcf


class KCFtoBedGraph:
    def __init__(self, args):
        self.input_file = args.input_file
        self.output_prefix = args.output_prefix
        self.sample_name = args.sample_name

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
