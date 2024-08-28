#!/usr/bin/env python
# Description: Convert IBSpy tsv file to kcf file

import os
from collections import defaultdict
from kcftools.utils.readers import read_tsv
from kcftools.utils.writers import write_kcf


class TSVtoKCF:
    def __init__(self, args):
        self.input_file = args.input
        self.output_file = args.output
        if args.sample:
            self.sample_name = args.sample
        else:
            self.sample_name = None

    def run(self):
        """
        Convert IBSpy tsv file to kcf file
        """
        if self.sample_name is None:
            sample_name = os.path.basename(self.input_file).split('.')[0]
        windows = defaultdict()
        with open(self.output_file, 'w') as o:
            for window in read_tsv(self.input_file, sample_name):
                windows[(window.chr_name, window.start, window.end)] = window
        write_kcf(self.output_file, windows, sample_name)
