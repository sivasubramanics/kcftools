import logging
from collections import defaultdict

from kcftools.data.IBS import IBS
from kcftools.utils.readers import read_kcf
from kcftools.utils.helper_functions import ibs_to_binary


class ExtractIBS:
    def __init__(self, args):
        self.input_file = args.input
        self.output_prefix = args.output
        self.sample_name = args.sample
        self.heatmap = args.heatmap

    def run(self):
        """
        Extract IBS regions from kcf file
        """
        windows, samples, misc_lines = read_kcf(self.input_file)
        logging.info(f'Extracting windows from {self.input_file}')
        if self.sample_name is not None:
            if type(self.sample_name) == str:
                samples = [self.sample_name]
            elif type(self.sample_name) == list:
                samples = self.sample_name

        for sample in samples:
            output_file = f'{self.output_prefix}.{sample}.tsv'
            fo = open(output_file, 'w')
            fo.write(f'id\tseqname\tstart\tend\tlength\ttotal_blocks\tibs_blocks\tmean_score\n')
            out_bed = f'{self.output_prefix}.{sample}.bed'
            fo_bed = open(out_bed, 'w')

            blocks = defaultdict()
            last_num = None
            for (seqname, start, end) in windows:
                block_num = windows[(seqname, start, end)].data[sample].ibs
                if block_num == 'N':
                    is_ibs = False
                    if last_num is not None:
                        blocks[last_num].add_window(start, end, windows[(seqname, start, end)].data[sample].score, is_ibs)
                else:
                    is_ibs = True
                    if block_num not in blocks:
                        blocks[block_num] = IBS(seqname)
                    blocks[block_num].add_window(start, end, windows[(seqname, start, end)].data[sample].score, is_ibs)
                    last_num = block_num

            for block_num in blocks:
                fo.write(f'{block_num}\t{blocks[block_num]}\n')
                for i in range(len(blocks[block_num].starts)):
                    fo_bed.write(f'{blocks[block_num].seqname}\t{blocks[block_num].starts[i]}\t{blocks[block_num].ends[i]}\t0\t+\n')

        if self.heatmap:
            prev_chrom = None
            num_window = 0
            for key in windows:
                chrom = key[0]
                num_window += 1
                if chrom != prev_chrom:
                    if prev_chrom is not None:
                        fo.close()
                    num_window = 1
                    fo = open(f'{self.output_prefix}.{chrom}.heatmap.tsv', 'w')
                    samples_line = '\t'.join(samples)
                    fo.write(f'window\t{samples_line}\n')
                ibs = "\t".join(map(ibs_to_binary, windows[key].get_value('IB', samples)))
                fo.write(f'{num_window}\t{ibs}\n')
                prev_chrom = chrom
            fo.close()
