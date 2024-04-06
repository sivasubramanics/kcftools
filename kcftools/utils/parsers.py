#!/usr/bin/env python3

from collections import defaultdict
from kcftools.utils.readers import read_kcf
from kcftools.utils.writers import write_kcf
from kcftools.data.Window import Window
from kcftools.utils.helper_functions import get_chromosomes, get_window


def recalculate(in_kcf, out_kcf):
    """
    Recalculate the scores in kcf file
    """
    windows, samples, misc_lines = read_kcf(in_kcf)
    write_kcf(out_kcf, windows, samples, misc_lines)


def convert_windows(windows, chr_lengths, o_win_size, n_win_size, kmer_size):
    """
    Convert window size
    """
    # get all the chromosomes from the windows dictionary
    chromosomes = get_chromosomes(windows)
    new_windows = defaultdict(dict)
    for seqname in chr_lengths:
        # check if the chromosome is in the fai file
        if seqname not in chromosomes:
            continue
        n_win_start = 0
        win_starts, win_ends = get_window(0, chr_lengths[seqname], n_win_size, kmer_size)
        for win_start, win_end in zip(win_starts, win_ends):
            new_window = Window(seqname, win_start, win_end)
            new_key = (str(seqname), str(win_start), str(win_end))
            o_win_starts, o_win_ends = get_window(n_win_start, win_end, o_win_size, kmer_size)
            for o_win_start, o_win_end in zip(o_win_starts, o_win_ends):
                key = (str(seqname), o_win_start, o_win_end)
                if key not in windows:
                    n_win_start = o_win_start
                    continue
                old_window = windows[key]
                samples = old_window.data.keys()
                new_window.add_total_kmers(old_window.total_kmers)
                for sample in samples:
                    new_window.add_data(sample, old_window.data[sample])
            new_windows[new_key] = new_window
    return new_windows


def update_info(line, sample_indices=None):
    """
    Update info field, it takes list of KCF fields and list of samples indices
    """
    attr_indices = []
    va = []
    ob = []
    di = []
    sc = []
    if sample_indices is None:
        sample_indices = range(6, len(line))
    for i in line[5].split(':'):
        if i in ['IB', 'VA', 'OB', 'DI', 'SC']:
            attr_indices.append(line[5].split(':').index(i))
    for i in sample_indices:
        va.append(float(line[i].split(':')[attr_indices[1]]))
        ob.append(float(line[i].split(':')[attr_indices[2]]))
        di.append(float(line[i].split(':')[attr_indices[3]]))
        sc.append(float(line[i].split(':')[attr_indices[4]]))
    return (f"MIN_SCORE={min(sc)};"
            f"MAX_SCORE={max(sc)};"
            f"MEAN_SCORE={round(sum(sc) / len(sc), 2)};"
            f"MIN_VA={min(va)};"
            f"MAX_VA={max(va)};"
            f"MEAN_VA={round(sum(va) / len(va), 2)};"
            f"MIN_OB={min(ob)};"
            f"MAX_OB={max(ob)};"
            f"MEAN_OB={round(sum(ob) / len(ob), 2)};")
