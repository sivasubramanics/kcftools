#!/usr/bin/env python3
# Description: This script is used to run IBSpy

import logging
import sys
import os

from kcftools.utils.helper_functions import run_cmd


class RunIBSpy:
    """
    Run IBSpy
    """
    def __init__(self, args):
        self.input = args.input
        self.reference = args.reference
        self.name = args.name
        self.threads = args.threads
        self.kmer = args.kmer
        self.window = args.window

    def run(self):
        """
        Run IBSpy
        """

        # check if kmc db is present
        kmer_db = f"{self.input}/kmc{self.kmer}"
        if not os.path.exists(kmer_db):
            logging.error(f"Kmer database not found: {kmer_db}")
            logging.error("Please run get_kmers command first")
            sys.exit(1)

        # check if reference fasta fai is present
        ref_fai = f"{self.reference}.fai"
        if not os.path.exists(ref_fai):
            logging.info(f"Creating reference fasta index: {ref_fai}")
            run_cmd(f"samtools faidx {self.reference}")

        ibspy_cmd = [f"IBScpp"]
        ibspy_cmd.append(f"-d {kmer_db}")
        ibspy_cmd.append(f"-r {self.reference}")
        ibspy_cmd.append(f"-p {self.threads}")
        ibspy_cmd.append(f"-k {self.kmer}")
        ibspy_cmd.append(f"-w {self.window}")
        ibspy_cmd.append(f"> {self.input}/{self.name}.variations.tsv")

        logging.info(f"Running IBSpy")
        run_cmd(ibspy_cmd)

