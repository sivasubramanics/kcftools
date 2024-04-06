#!/usr/bin/env python3

"""
Increase window size
"""

import logging
from kcftools.utils.readers import read_kcf, get_seq_lengths, get_current_window_size
from kcftools.utils.writers import write_kcf



class IncreaseWindow:
    """
    Increase window size
    """
    def __init__(self, args):
        self.input_file = args.input_file
        self.output_file = args.output_file
        self.window_size = args.window_size
        self.kmer_size = args.kmer_size
        self.length_file = args.length_file

    def run(self):
        """
        Increase window size
        """
        windows, samples, misc_lines = read_kcf(self.input_file)
        seq_lengths = get_seq_lengths(self.length_file)
        current_window_size = get_current_window_size(windows, seq_lengths)
        logging.info(f'Current window size: {current_window_size}')
        windows = convert_windows(windows, seq_lengths, current_window_size, int(window_size), int(kmer_size))
        logging.info(f'New window size: {window_size}')
        write_kcf(output_file, windows, samples, misc_lines)



