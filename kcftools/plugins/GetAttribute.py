import logging
import sys


class GetAttribute:
    def __init__(self, args):
        self.input = args.input
        self.output = args.output
        self.samples = args.samples
        self.attribute = args.attribute

    def run(self):
        """
        Extract scores from kcf file for a given list of samples
        """
        if self.attribute not in ['VA', 'OB', 'DI', 'SC']:
            logging.info(f'Error: {self.attribute} is not a valid attribute')
            sys.exit(1)
        sample_indices = []
        win = 0
        if type(self.samples) == str:
            samples = [self.samples]
        elif type(self.samples) == list:
            samples = self.samples
        elif self.samples is None:
            samples = []
        else:
            logging.error('Error: samples should be a string or a list of strings')
            sys.exit(1)
        fi = open(self.input, 'r')
        fo = open(self.output, 'w')
        fm = open(self.output + '.map', 'w')
        logging.info(f'Extracting {self.attribute} from {self.input}')
        for line in fi:
            if line.startswith('##'):
                continue
            line = line.strip().split('\t')
            if line[0].startswith('#CHROM'):
                if not samples:
                    samples = line[6:]
                for sample in samples:
                    if sample in line:
                        sample_indices.append(line.index(sample))
                fm.write(f'WIN\t{line[0]}\t{line[1]}\t{line[2]}\n')
                fo.write('WIN')
                for i in sample_indices:
                    fo.write(f'\t{line[i]}')
                fo.write('\n')
                continue
            win += 1
            attribute_idx = line[5].split(':').index(self.attribute)
            fm.write(f'{win}\t{line[0]}\t{line[1]}\t{line[2]}\n')
            fo.write(f'{win}')
            for i in sample_indices:
                fo.write(f'\t{line[i].split(":")[attribute_idx]}')
            fo.write('\n')
        fi.close()
        fo.close()
