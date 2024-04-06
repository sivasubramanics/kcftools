#!/usr/bin/env python3
# Description: Contains the Window class

from collections import defaultdict


class Window:
    """
    Window class to store the values for a window
    """

    def __init__(self, chr_name, start, end):
        self.chr_name = chr_name
        self.start = start
        self.end = end
        self.total_kmers = 0
        self.data = defaultdict()
        self.blocks = defaultdict()

    def set_total_kmers(self, total_kmers):
        self.total_kmers = total_kmers

    def add_total_kmers(self, total_kmers):
        self.total_kmers += total_kmers

    def set_value(self, sample_names, attribute, value):
        if type(sample_names) == str:
            sample_names = [sample_names]
        for sample_name in sample_names:
            if sample_name not in self.data:
                self.data[sample_name] = Data(0, 0, 0, self.total_kmers)
            if attribute == 'VA':
                self.data[sample_name].va = value
            elif attribute == 'OB':
                self.data[sample_name].ob = value
            elif attribute == 'DI':
                self.data[sample_name].di = value

    def set_data(self, sample_name, data):
        if sample_name not in self.data:
            self.data[sample_name] = Data('N', 0, 0, 0, self.total_kmers)
        self.data[sample_name] = data

    def add_data(self, sample_name, data):
        if sample_name not in self.data:
            self.data[sample_name] = Data('N', 0, 0, 0, 0)
        self.data[sample_name].tot += data.tot
        self.data[sample_name].ob += data.ob
        self.data[sample_name].va += data.va
        self.data[sample_name].di += data.di

    def get_value(self, attribute, sample_names=None):
        if sample_names is None:
            sample_names = self.data.keys()
        if type(sample_names) == str:
            sample_names = [sample_names]
        values = []
        for sample_name in sample_names:
            if attribute == 'VA':
                values.append(self.data[sample_name].va)
            elif attribute == 'OB':
                values.append(self.data[sample_name].ob)
            elif attribute == 'DI':
                values.append(self.data[sample_name].di)
            elif attribute == 'SC':
                values.append(self.data[sample_name].score)
            elif attribute == 'IB':
                values.append(self.data[sample_name].ibs)
            else:
                sys.exit(f'Error: Attribute {attribute} is not valid')
        return values

    def get_window(self, sample_names=None):
        if sample_names is None:
            sample_names = self.data.keys()
        vcf_list = []
        if type(sample_names) == str:
            sample_names = [sample_names]
        for sample_name in sample_names:
            vcf_list.append(f"{self.data[sample_name]}")
        return '\t'.join(vcf_list)

    def get_info(self, sample_names=None):
        if sample_names is None:
            sample_names = self.data.keys()
        if type(sample_names) == str:
            sample_names = [sample_names]
        scores = self.get_value('SC', sample_names)
        va = self.get_value('VA', sample_names)
        ob = self.get_value('OB', sample_names)
        di = self.get_value('DI', sample_names)
        return (f"MIN_SCORE={min(scores)};"
                f"MAX_SCORE={max(scores)};"
                f"MEAN_SCORE={round(sum(scores) / len(scores), 2)};"
                f"MIN_VA={min(va)};"
                f"MAX_VA={max(va)};"
                f"MEAN_VA={round(sum(va) / len(va), 2)};"
                f"MIN_OB={min(ob)};"
                f"MAX_OB={max(ob)};"
                f"MEAN_OB={round(sum(ob) / len(ob), 2)};")

    def __str__(self):
        return f'{self.chr_name}\t{self.start}\t{self.end}'


class Data:
    """
    Data class to store the values for a sample
    """

    def __init__(self, ibs, ob, va, di, tot):
        self.ob = ob
        self.va = va
        self.di = di
        self.tot = tot
        self.ibs = ibs

    @property
    def score(self):
        if self.tot == 0:
            return 0.0000
        if self.ob == 0:
            return 0.0000
        _score = round((self.ob - self.va) / self.tot, 4)
        if _score < 0.0001:
            return 0.0000
        return _score

    def __str__(self):
        return f'{self.ibs}:{self.va}:{self.ob}:{self.di}:{self.score}'
