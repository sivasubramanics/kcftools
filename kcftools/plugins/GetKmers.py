#!/usr/bin/env python3
# Description: This script is used to run IBSpy

import logging
import sys
import os
import time
from typing import List

from kcftools.utils.helper_functions import run_cmd
from kcftools.utils.helper_functions import is_installed
from kcftools._defaults import ibspy_tools


class GetKmers:
    def __init__(self, args):
        print(args)
        if args.input:
            self.input = args.input
        self.output = args.output
        if args.input_list:
            self.input_list = args.input_list
        else:
            self.input_list = None
        if not self.input_list and not self.input:
            logging.error("Please provide input or input_list")
            exit(1)
        self.threads = args.threads
        self.kmer = args.kmer
        self.format = args.format
        self.mem = args.mem
        self.timestamp = str(time.time()).replace(".", "")

    def run(self):
        """
        Run IBSpy
        """
        # create random variable using timestamp to create temporary file

        for tool in ibspy_tools:
            if not is_installed(tool):
                logging.error(f"{tool} is not installed. Please install it and try again")
                sys.exit(1)

        # create output directory if not exists
        if not os.path.exists(self.output):
            logging.info(f"Creating output directory: {self.output}")
            os.makedirs(self.output, exist_ok=True)
            os.makedirs(f"{self.output}/tmp", exist_ok=True)

        # run kmc
        self.run_kmc()

        # run kmers_add_strand_information
        self.run_kmers_add_strand_information()

    def run_kmc(self):
        """
        Run kmc for cannonized and all kmers
        """
        kmc_cmd = [f"kmc"]
        # threads
        kmc_cmd.append(f"-t{self.threads}")
        # memory
        kmc_cmd.append(f"-m{self.mem}")
        # kmersize
        kmc_cmd.append(f"-k{self.kmer}")
        # minimum count
        kmc_cmd.append("-ci0")
        # hide progress
        kmc_cmd.append("-hp")
        # input format
        if self.format == 'fa':
            kmc_cmd.append("-fa")
        if self.format == 'fq':
            kmc_cmd.append("-fq")
        if self.format == 'fm':
            kmc_cmd.append("-fm")

        if self.input_list:
            kmc_ifiles = f"@{self.input_list}"
        if self.input:
            if len(self.input) == 1:
                kmc_ifiles = self.input[0]
            else:
                # If multiple input files are provided, create temporary file inside output directory with list of input files
                # make the filename starts with timestamp to avoid conflicts
                with open(f"{self.output}/{self.timestamp}_input_list.txt", 'w') as file:
                    for i in self.input:
                        # absolute path of i
                        ipath = os.path.abspath(i)
                        file.write(f"{ipath}\n")
                kmc_ifiles = f"@{self.output}/{self.timestamp}_input_list.txt"

        # input files
        kmc_cmd.append(kmc_ifiles)

        # output file
        kmc_cmd.append(f"{self.output}/canon")

        # temporary directory
        kmc_cmd.append(f"{self.output}/tmp")

        # run kmc
        stdout, stderr = run_cmd(kmc_cmd, "Running kmc for cannonized kmers")
        # write stdout to log file
        with open(f"{self.output}/canon.log", 'w') as log:
            log.write(stdout)

        # remove last 3 from kmc_cmd to run kmc for all kmers
        print(kmc_cmd)
        kmc_cmd = kmc_cmd[:-3]
        kmc_cmd.append("-b")
        print(kmc_cmd)
        kmc_cmd.append(f"{kmc_ifiles}")
        kmc_cmd.append(f"{self.output}/all")
        kmc_cmd.append(f"{self.output}/tmp")

        # run kmc for all kmers
        stdout, stderr = run_cmd(kmc_cmd, "Running kmc for all kmers")
        # write stdout to log file
        with open(f"{self.output}/all.log", 'w') as log:
            log.write(stdout)

    def run_kmers_add_strand_information(self):
        """
        Run kmers_add_strand_information
        """
        canon_kmers = f"{self.output}/canon"
        all_kmers = f"{self.output}/all"

        kmc_add_strand_cmd = [f"kmers_add_strand_information"]
        kmc_add_strand_cmd.append(f"-c {canon_kmers}")
        kmc_add_strand_cmd.append(f"-n {all_kmers}")
        kmc_add_strand_cmd.append(f"-k {self.kmer}")
        kmc_add_strand_cmd.append(f"-o {self.output}/kmc{self.kmer}")

        # run kmers_add_strand_information
        stdout, stderr = run_cmd(kmc_add_strand_cmd, "Running kmers_add_strand_information (this may take a while)")
        # write stdout to log file
        with open(f"{self.output}/kmc{self.kmer}.log", 'w') as log:
            log.write(stdout)
        logging.info(f"Output is written to {self.output}/kmc{self.kmer}")




