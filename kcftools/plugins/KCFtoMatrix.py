from kcftools.utils.helper_functions import convert_tsv2RData, transpose_matrix
from kcftools.utils.readers import read_kcf
from kcftools.utils.writers import write_matrix




class KCFtoMatrix:
    def __init__(self, args):
        self.input_file = args.input_file
        self.output_prefix = args.output_prefix
        self.sample = args.sample
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
        write_matrix(windows, self.output_prefix, samples, self.allele_a_cutoff, self.allele_b_cutoff)
        # transpose_gt_matrix(f"{outprefix}.matrix.tr.tsv", f"{outprefix}.matrix.tsv")
        transpose_matrix(f"{self.output_prefix}.matrix.tr.tsv", f"{self.output_prefix}.matrix.tsv")
        convert_tsv2RData(f"{self.output_prefix}.matrix.tsv", f"{self.output_prefix}.map.tsv")


