

import sys

from kcftools.utils.helper_functions import list_to_str


class SplitKCF:
    def __init__(self, args):
        self.in_kcf = args.in_kcf
        self.output_prefix = args.output_prefix
        self.sample_name = args.sample_name

    def run(self, out_prefix=None):
        """
        Split kcf file by chromosome for given set of samples if list of samples provided
        """
        prev_chrom = None
        data_indices = [0, 1, 2, 3, 4, 5]
        misc_lines = []
        with open(self.in_kcf, 'r') as f:
            for line in f:
                if line.startswith('##'):
                    misc_lines.append(line.strip())
                    continue
                if line.startswith('#CHROM'):
                    misc_lines.append(f"##CMD: {' '.join(sys.argv)}")
                    in_samples = line.strip().split('\t')[6:]
                    if samples is None:
                        samples = in_samples
                    if type(samples) == str:
                        samples = [samples]
                    for sample in samples:
                        data_indices.append(in_samples.index(sample) + 6)
                        if sample not in in_samples:
                            sys.exit(f'Error: {sample} not found in {self.in_kcf}')
                    misc_lines.append(f'#CHROM\tSTART\tEND\tTOTAL_KMER\tINFO\tFORMAT\t{list_to_str(samples)}')
                    continue
                else:
                    line = line.strip().split('\t')
                    out_line = [line[i] for i in data_indices]
                    if chrs is not None:
                        if line[0] not in chrs:
                            continue
                    if prev_chrom is None or line[0] != prev_chrom:
                        print_log(f'Writing {out_prefix}.{line[0]}.kcf')
                        fo = open(f'{out_prefix}.{line[0]}.kcf', 'w')
                        fo.write(f'{list_to_str(misc_lines, NEWLINE)}\n')
                    fo.write(f'{list_to_str(out_line)}\n')
                    prev_chrom = line[0]
