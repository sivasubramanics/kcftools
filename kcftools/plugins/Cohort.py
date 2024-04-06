#!/usr/bin/env python3
# Description: Contains the Cohort class

from collections import defaultdict
from kcftools.utils.readers import read_kcf
from kcftools.utils.writers import write_kcf

class Cohort:
    """
    Combine multiple kcf files
    """
    def __init__(self, args):
        self.input_files = args.input_files
        self.output_file = args.output_file

    def run(self):
        """
        Combine multiple kcf files
        """
        windows = defaultdict()
        samples = []
        with open(self.input_files, 'r') as f:
            for input_file in f:
                input_file = input_file.strip()
                windows, samples, misc_lines = read_kcf(input_file, windows, samples)
        write_kcf(self.output_file, windows, samples, misc_lines)

