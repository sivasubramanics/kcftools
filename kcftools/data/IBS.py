

class IBS:
    def __init__(self, seqname):
        self.seqname = seqname
        self.starts = []
        self.ends = []
        self.scores = []
        self.is_ibs = []

    def add_window(self, start, end, score, is_ibs):
        self.starts.append(start)
        self.ends.append(end)
        self.scores.append(score)
        self.is_ibs.append(is_ibs)

    def remove_tail_na(self):
        # remove the tailing False from the list of is_ibs
        # count number of False from the end of the list
        tail_na_count = 0
        for is_ibs in reversed(self.is_ibs):
            if is_ibs:
                break
            tail_na_count += 1
        if tail_na_count == 0:
            return
        self.starts = self.starts[:-tail_na_count]
        self.ends = self.ends[:-tail_na_count]
        self.scores = self.scores[:-tail_na_count]
        self.is_ibs = self.is_ibs[:-tail_na_count]

    def __str__(self):
        self.remove_tail_na()
        return (f'{self.seqname}'
                f'\t{self.starts[0]}'
                f'\t{self.ends[-1]}'
                f'\t{self.ends[-1]-self.starts[0]}'
                f'\t{len(self.starts)}'
                f'\t{len([x for x in self.is_ibs if x])}'
                f'\t{round(sum(self.scores) / len(self.scores), 2)}')