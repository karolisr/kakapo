# -*- coding: utf-8 -*-

"""ORF Finder"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

from kakapo.seq import reverse_complement


def find_stop_codon_free_zones(seq, forward_frame, min_len, stop_codons):  # noqa
    assert forward_frame > 0 and forward_frame < 4

    seq = seq.upper()

    offset = forward_frame - 1
    # Adjust sequence so it is in the requested frame
    in_frame = seq[offset:len(seq)]
    # Split sequence into codons
    codons = [(in_frame[i:i + 3]) for i in range(0, len(in_frame), 3)]
    codons = [x for x in codons if len(x) == 3]
    # Identify stop codon positions in codons
    idx_stop = [i for i, x in enumerate(codons) if x in stop_codons]
    idx_stop = [-1] + idx_stop + [len(codons)]
    # Determine interval sizes between stop codons
    il = [(idx_stop[i + 1] - idx_stop[i]) for i in range(len(idx_stop) - 1)]
    il = [idx_stop[0]] + il
    # Indexes of reverse-sorted interval sizes between stop codons
    idx_int = []
    for l in sorted(il, reverse=True):
        idx = il.index(l)
        # index() finds the first occurence of a value in the list. Replacing
        # the found value prevents it from being considered again.
        il[idx] = ''
        idx_int.append(idx)

    # Stop codon-free zone coordinates sorted by their size from largest to
    # smallest
    zones = []
    for idx_zone_begin in idx_int:
        idx_zone_end = idx_zone_begin + 2
        zone = idx_stop[idx_zone_begin:idx_zone_end]
        if len(zone) != 2:
            continue
        # Convert back to nucleotide coordinates
        zone_nt = [(zone[0] * 3) + offset + 3, (zone[1] * 3) + offset]
        if zone_nt[1] - zone_nt[0] < min_len:
            continue
        zones.append(zone_nt)

    return zones


def find_start_codons(seq, start_codons):
    """Argument seq is assumed to be in frame."""
    assert type(start_codons) in (list, tuple, set)

    seq = seq.upper()

    # Split sequence into codons
    codons = [(seq[i:i + 3]) for i in range(0, len(seq), 3)]

    # Identify start codon positions in codons
    idx_start = [i for i, x in enumerate(codons) if x in start_codons]
    # Convert back to nucleotide coordinates
    idx_start_nt = [i * 3 for i in idx_start]

    return idx_start_nt


def find_orfs(seq, frame, min_len, stop_codons, start_codons=None,
              include_terminal_codon=True):  # noqa
    assert frame > -4 and frame != 0 and frame < 4
    assert start_codons is None or type(start_codons) in (list, tuple, set)

    seq = seq.upper()

    if frame < 0:
        seq = reverse_complement(seq)

    zones = find_stop_codon_free_zones(seq, abs(frame), min_len, stop_codons)

    if start_codons is not None:
        for zone in zones:
            seq_in_zone = seq[zone[0]:zone[1]]
            idx_start_nt = find_start_codons(seq_in_zone, start_codons)
            if len(idx_start_nt) > 0:
                # Considers only the earliest start codon: idx_start_nt[0]
                # Effectively this returns the longest ORF for the zone
                zone[0] = zone[0] + idx_start_nt[0]

    # Include stop codon?
    if include_terminal_codon is True:
        zones_with_terminal_codon = []
        for zone in zones:
            beg = zone[0]
            end = zone[1]
            end_with_term = end + 3
            if len(seq) >= end_with_term:
                end = end_with_term
            zones_with_terminal_codon.append([beg, end])

        zones = zones_with_terminal_codon

    if frame < 0:
        zones = [[len(seq) - x[1], len(seq) - x[0]] for x in zones]

    return zones


def find_longest_orfs(seq, frame, min_len, stop_codons, start_codons=None,
                      include_terminal_codon=True, ret_max=1):  # noqa

    orfs = find_orfs(seq, frame, min_len, stop_codons, start_codons,
                     include_terminal_codon)

    if len(orfs) == 0:
        return None

    lengths = [x[1] - x[0] for x in orfs]

    # Indexes of reverse-sorted interval sizes between stop codons
    idx_len = []
    for l in sorted(lengths, reverse=True):
        idx = lengths.index(l)
        # index() finds the first occurence of a value in the list. Replacing
        # the found value prevents it from being considered again.
        lengths[idx] = ''
        idx_len.append(idx)

    orfs_sorted = [orfs[o] for o in idx_len]

    return orfs_sorted[0:ret_max]


def find_orf_for_blast_hit(seq, frame, hit_start, hit_end, hit_length_adjust,
                           stop_codons, start_codons=None,
                           include_terminal_codon=True):  # noqa

    # If BLAST hit start occurs before the first start codon and the ORF
    # without start codon starts before the BLAST hit, return the ORF without
    # start codon.

    # min_len = hit_end - hit_start - int(hit_reduce)
    hit_len = hit_end - hit_start
    min_len = hit_len * hit_length_adjust

    if start_codons is None:
        sc_to_consider = [None]
    else:
        sc_to_consider = [['ATG']]
        start_codons = list(set(start_codons).difference(['ATG']))
        if len(start_codons) == 0:
            sc_to_consider.append(None)
        else:
            sc_to_consider.append(start_codons)
            sc_to_consider.append(None)

    orf = None
    for sc in sc_to_consider:
        orfs = find_longest_orfs(seq, frame, min_len, stop_codons, sc,
                                 include_terminal_codon=include_terminal_codon,
                                 ret_max=100)

        if orfs is None:
            continue

        for orf_to_consider in orfs:

            orf_start = orf_to_consider[0]
            orf_end = orf_to_consider[1]

            adj = (hit_len - min_len) / 2

            if orf_start <= (hit_start + adj) and orf_end >= (hit_end - adj):
                orf = orf_to_consider
                break

        if orf is not None:
            break

    return orf
