#!/usr/bin/env python3

"""
Run engine for IBS analysis
"""

import logging
import sys
import time


def run(args):
    """
    Run IBSpy
    """
    start_time = time.time()
    args = ArgumentParser(args)
    if args.command == 'get_kmers':
        from kcftools.plugins.GetKmers import GetKmers
        GetKmers.run(args)
    elif args.command == 'run_ibspy':
        from kcftools.plugins.RunIBSpy import RunIBSpy
        RunIBSpy.run(args)
    elif args.command == 'tsv2kcf':
        from kcftools.plugins.TSVtoKCF import TSVtoKCF
        TSVtoKCF.run(args)
    elif args.command == 'cohort':
        from kcftools.plugins.Cohort import Cohort
        Cohort.run(args)
    elif args.command == 'increase_window':
        from kcftools.plugins.IncreaseWindow import IncreaseWindow
        IncreaseWindow.run(args)
    elif args.command == 'find_IBS':
        from kcftools.plugins.FindIBS import FindIBS
        FindIBS.run(args)
    elif args.command == 'extract':
        from kcftools.plugins.ExtractIBS import ExtractIBS
        ExtractIBS.run(args)
    elif args.command == 'kcf2bedgraph':
        from kcftools.plugins.KCFtoBedGraph import KCFtoBedGraph
        KCFtoBedGraph.run(args)
    elif args.command == 'kcf2matrix':
        from kcftools.plugins.KCFtoMatrix import KCFtoMatrix
        KCFtoMatrix.run(args)
    elif args.command == 'split_kcf':
        from kcftools.plugins.SplitKCF import SplitKCF
        SplitKCF.run(args)
    elif args.command == 'concat':
        from kcftools.plugins.ConcatKCF import ConcatKCF
        ConcatKCF.run(args)
    elif args.command == 'extract_samples':
        from kcftools.plugins.ExtractSamples import ExtractSamples
        ExtractSamples.run(args)
    elif args.command == 'get_attr':
        from kcftools.plugins.GetAttribute import GetAttribute
        GetAttribute.run(args)
    elif args.command == 'extract_score_matrix':
        from kcftools.plugins.IBStoMatrix import IBStoMatrix
        IBStoMatrix.run(args)
    elif args.command == 'count_genes':
        from kcftools.plugins.CountGenes import CountGenes
        CountGenes.run(args)
    elif args.command == 'recalculate':
        from kcftools.plugins.Recalculate import Recalculate
        Recalculate.run(args)
    else:
        logging.error(f"Unknown command: {args.command}")
        sys.exit(1)

    logging.info("Done")
    end_time = time.time()
    runtime = end_time - start_time
    logging.info(f"Runtime: {runtime:.2f} seconds")

