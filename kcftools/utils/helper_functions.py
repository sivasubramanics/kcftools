#!/usr/bin/env python3
# Description: contains helper functions used in the project

import logging
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
import time
import pandas as pd


def run_cmd(cmd, message=None):
    """
    Run a command and return stdout, stderr as string
    """
    try:
        if message:
            logging.info(message)
        # if cmd is list, convert it to string
        if isinstance(cmd, list):
            cmd = ' '.join(cmd)
        logging.info(f"CMD: {cmd}")
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        stdout = stdout.decode('utf-8').strip()
        stderr = stderr.decode('utf-8').strip()
        run_code = process.returncode
        if run_code != 0:
            logging.error(f"Error running the command: {cmd}")
            logging.error(f"Error: {stderr}")
            sys.exit(1)
        return stdout, stderr
    except Exception as e:
        logging.error(f"Error running the command: {cmd}")
        logging.error(f"Error: {e}")
        sys.exit(1)


def is_installed(tool):
    """
    Check if the tool is accessible in the system. Return True if installed, False otherwise
    """
    # cmd = f"which {tool}"
    # stdout, stderr = run_cmd(cmd)
    # if stdout:
    #     return True
    # return False
    if shutil.which(tool):
        return True
    return False


def duplicates(it):
    """
    Check if there are duplicates in a list and return the duplicates
    """
    return [item for item, count in Counter(it).items() if count > 1]


def get_chromosomes(windows):
    """
    Return a list of chromosomes in the windows
    """
    chromosomes = []
    for window in windows:
        if window[0] not in chromosomes:
            chromosomes.append(window[0])
    return chromosomes


def get_current_window_size(windows, chr_lengths):
    """
    Get current window size
    """
    window_size = 0
    for (seqname, start, end) in windows:
        if window_size == 0:
            window_size = end - start
        if end == chr_lengths[seqname]:
            continue
        if window_size != end - start:
            sys.exit(f'Error: Window size is not consistent in input file.')
    return int(window_size)


def get_window(start, end, window_size, kmer_size):
    """
    return a list of windows
    """
    win_start = []
    win_end = []
    t_win_start = start
    t_win_end = start + window_size
    while t_win_end < end:
        win_start.append(t_win_start)
        win_end.append(t_win_end)
        t_win_start = t_win_end - kmer_size
        t_win_end = t_win_start + window_size
    if t_win_end >= end:
        win_start.append(t_win_start)
        win_end.append(end)
    return win_start, win_end


def get_seq_lengths(length_file):
    """
    Read sequence length file and return a dictionary
    """
    chr_lengths = defaultdict()
    with open(length_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            chr_lengths[line[0]] = int(line[1])
    return chr_lengths


def ibs_to_binary(in_str):
    """
    Convert IBS to binary
    """
    if in_str == 'N':
        return '0'
    else:
        return '1'


def monomorphic(genotypes):
    """
    Check if the genotypes are monomorphic
    """
    # delete the Ns in geneotypes
    genotypes = [x for x in genotypes if x != 'N']
    genotypes = list(set(genotypes))
    if len(genotypes) == 1:
        return True
    else:
        return False


def convert_tsv2RData(*tsv_list):
    """
    Convert tsv file to RData file
    """
    # if Rscript is not installed, exit
    if not shutil.which('Rscript'):
        sys.exit('Error: Rscript is not installed')

    for in_tsv in tsv_list:
        out_rdata = f"{in_tsv}.RData"
        logging.info(f'Converting {in_tsv} to {out_rdata}')
        rscript = f"{os.path.dirname(os.path.realpath(__file__))}/{time.strftime('%Y%m%d-%H%M%S')}.tsv2RData.R"
        rscript_fo = open(rscript, 'w')
        if in_tsv.endswith('.matrix.tsv'):
            rscript_fo.write(f"myGD <- read.table(\"{os.path.abspath(in_tsv)}\", header=TRUE)\n")
            rscript_fo.write(f"save(myGD, file=\"{os.path.abspath(out_rdata)}\")\n")
        elif in_tsv.endswith('.map.tsv'):
            rscript_fo.write(f"myGM <- read.table(\"{os.path.abspath(in_tsv)}\", header=TRUE)\n")
            rscript_fo.write(f"save(myGM, file=\"{os.path.abspath(out_rdata)}\")\n")
        rscript_fo.close()
        cmd = f"Rscript {rscript}"
        os.system(cmd)
        os.remove(rscript)


def transpose_matrix(in_file, out_file):
    """
    Transpose matrix
    """
    logging.info(f'Transposing {in_file}')
    df = pd.read_csv(in_file, sep='\t')
    df = df.T
    df.to_csv(out_file, sep='\t', header=False)


def list_to_str(in_list, sep='\t'):
    """
    Convert list to string
    """
    return sep.join(map(str, in_list))
