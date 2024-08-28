import logging

from kcftools.utils.helper_functions import convert_tsv2RData, transpose_matrix
from kcftools.utils.readers import read_kcf
from kcftools.utils.writers import write_matrix


class KCFtoMatrix:
    """
    Convert kcf file to genotype matrix and genotype map file.
    """

    def __init__(self, args):
        self.input_file = args.input
        self.output_prefix = args.output
        if args.sample:
            self.sample = args.sample
        else:
            self.sample = None
        self.allele_a_cutoff = args.allele_a_cutoff
        self.allele_b_cutoff = args.allele_b_cutoff

    def run(self):
        """
        Convert kcf file to genotype matrix and genotype map file.
        the map file will have the window name, chromosome and the position in TSV format
        """
        windows, samples, misc_lines = read_kcf(self.input_file)
        if self.sample is not None:
            if type(self.sample) == str:
                samples = [self.sample]
            elif type(self.sample) == list:
                samples = self.sample
        if len(samples) == 1:
            logging.error("Only one sample found. Please provide more than one sample")
            return
        write_matrix(windows, self.output_prefix, samples, self.allele_a_cutoff, self.allele_b_cutoff)
        # TODO: as of now i am using a crude way of transposing the matrix using Pandas. For larger files, this may not be efficient.
        # TODO: Need to find a better way to transpose the matrix. [not on priority as of now]
        # transpose_gt_matrix(f"{outprefix}.matrix.tr.tsv", f"{outprefix}.matrix.tsv")
        transpose_matrix(f"{self.output_prefix}.matrix.tr.tsv", f"{self.output_prefix}.matrix.tsv")
        convert_tsv2RData(f"{self.output_prefix}.matrix.tsv", f"{self.output_prefix}.map.tsv")


