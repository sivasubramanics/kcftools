#!/usr/bin/env python3

from kcftools.utils.parsers import recalculate


class Recalculate:
    def __init__(self, args):
        self.input = args.input
        self.output = args.output

    def run(self):
        """
        Run Recalculate attributes in KCF
        """
        recalculate(self.input, self.output)