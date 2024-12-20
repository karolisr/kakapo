"""ORF related code."""

from functools import partial, reduce
from itertools import accumulate, chain, compress, dropwhile, groupby, starmap
from operator import add, contains, itemgetter, mul, ne, not_, sub
from statistics import stdev

from kakapo.tools.seq import reverse_complement
from kakapo.utils.misc import overlap

CONTEXT_NT_CODES = str.maketrans('ACGTUBDHKMNRSVWY', '0123344444444444')


def get_codons(seq):
    """Split seq into codons."""
    reduce_add = partial(reduce, add)
    return tuple(map(reduce_add, zip(*[iter(seq)] * 3)))


def separate_stop_codons(codons, stop_codons):
    test = partial(contains, stop_codons)
    return tuple(tuple(v) for k, v in groupby(codons, test))


def keep_orf_codons(sep_codons, start_codons):
    test = partial(lambda *x: not_(contains(*x)), start_codons)
    return tuple(map(lambda x: tuple(dropwhile(test, x)), sep_codons))


def get_orf_coords(seq, start_codons, stop_codons):
    sep_codons = separate_stop_codons(get_codons(seq), stop_codons)

    orf_len = tuple(map(len, keep_orf_codons(sep_codons, start_codons)))
    orf_idx = tuple(map(partial(ne, 0), orf_len))

    orf_end = tuple(compress(accumulate(map(len, sep_codons)), orf_idx))
    orf_beg = starmap(sub, zip(orf_end, compress(orf_len, orf_idx)))

    return tuple(zip(map(partial(mul, 3), orf_beg),
                     map(partial(mul, 3), orf_end)))


def get_orf_coords_for_forward_frame(seq, start_codons, stop_codons,
                                     min_len=100, frame=1):
    assert type(frame) is int
    assert frame >= 1
    if len(seq) < min_len:
        return tuple()
    offset = frame - 1
    orf_coords = get_orf_coords(seq[offset:], start_codons, stop_codons)
    adjust = partial(add, offset)
    coords = tuple(map(lambda x: tuple(map(adjust, x)), orf_coords))
    coords_filtered = tuple(filter(lambda x: x[1] - x[0] >= min_len, coords))
    if len(coords_filtered) > 0:
        used_start_codons = tuple(map(lambda x: x[0] // 3, coords))
        seq = seq[:offset % 3] + reduce(add, ['...' if i in used_start_codons else x
                                              for i, x in enumerate(get_codons(seq[offset % 3:]))])
        new_start = min(tuple(chain.from_iterable(coords_filtered))) + 4
        coords_filtered += get_orf_coords_for_forward_frame(
            seq, start_codons, stop_codons, min_len, frame=new_start)
    return coords_filtered


def get_orf_coords_for_frames(seq, start_codons, stop_codons, min_len=100,
                              frames=(-3, -2, -1, 1, 2, 3)):

    valid_frames = {-3, -2, -1, 1, 2, 3}
    assert valid_frames | set(frames) == valid_frames

    sl = len(seq)
    results = []

    for frame in frames:
        s = seq
        if frame < 0:
            s = reverse_complement(seq)
        coords = get_orf_coords_for_forward_frame(s, start_codons, stop_codons,
                                                  min_len, abs(frame))
        if frame < 0:
            coords = tuple(map(lambda x: (sl - x[1], sl - x[0]), coords))

        result = dict()
        result['frame'] = frame
        result['orf_coords'] = coords

        results.append(result)

    return results


def geometric_mean(x):
    x = tuple(x)
    y = reduce(mul, map(lambda i: 1 + i, x), 1)
    return (y ** (1 / len(x))) - 1


def start_codon_score_partial(seq, context):
    size = min(len(seq), len(context))
    if size == 0:
        return 0
    stdev_gm = partial(stdev, xbar=geometric_mean(map(max, context)))
    weights = tuple(map(stdev_gm, context))
    seq_encoded = seq.translate(CONTEXT_NT_CODES)
    context = [list(x) + [0] for x in context]
    data = tuple(zip(map(int, seq_encoded), context, weights))
    max_score = geometric_mean(map(lambda x: max(x[1]) * x[2], data))
    score = geometric_mean(starmap(lambda x, y, z: y[x] * z, data))
    return score / max_score


def start_codon_score(seq, idx, context_l, context_r):
    seq_l = seq[0:idx:][::-1]
    seq_r = seq[idx + 3:]
    q_l = start_codon_score_partial(seq_l, context_l)
    q_r = start_codon_score_partial(seq_r, context_r)
    if q_l == 0 or q_r == 0:
        return max(q_l, q_r)
    return (q_l * q_r) ** 0.5


def find_orf_for_blast_hit(seq, frame, hit_start, hit_end, start_codons,
                           stop_codons, context_l, context_r, min_overlap,
                           min_len, max_len, allow_no_strt_cod,
                           allow_no_stop_cod):

    seq = str(seq)

    assert type(frame) is int

    log_str = ''

    seq_len = len(seq)
    seq_local = seq
    min_len_to_consider = 200

    # FixMe: very long sequence (genomic) hack
    padding = 10000
    max_hit_len = padding
    seq_temp = seq
    seq_trim_start = 0
    seq_trim_end = seq_len
    ###

    if frame > 0:
        hit_len = hit_end - hit_start
        frames = [1, 2, 3]

    else:
        seq_local = reverse_complement(seq)
        hit_start_temp = hit_start
        hit_start = hit_end
        hit_end = hit_start_temp
        hit_len = hit_end - hit_start
        frames = [-1, -2, -3]

    # FixMe: very long sequence (genomic) hack
    if hit_start - padding >= 0:
        seq_trim_start = hit_start - padding
    if hit_end + padding <= seq_len:
        seq_trim_end = hit_end + padding
    seq_temp = seq[seq_trim_start:seq_trim_end]
    ###

    # print('len(seq)', len(seq))
    # print('len(seq_temp)', len(seq_temp))
    # print('frame', frame)
    # print(hit_start, hit_end)

    # FixMe: very long sequence (genomic) hack
    results = list()
    if hit_len < max_hit_len:
        results = get_orf_coords_for_frames(seq_temp, start_codons, stop_codons,
                                            min_len_to_consider, frames=frames)
    ###

    orfs = list()

    for r in results:

        orf_frame = r['frame']
        orf_coords = r['orf_coords']

        for loc in orf_coords:

            # FixMe: very long sequence (genomic) hack
            loc = list(loc)
            loc[0] = loc[0] + seq_trim_start
            loc[1] = loc[1] + seq_trim_start
            ###

            orf_beg = loc[0]
            orf_end = loc[1]

            if frame < 0:
                orf_beg = seq_len - loc[1]
                orf_end = seq_len - loc[0]

            sc_score = start_codon_score(seq_local, orf_beg,
                                         context_l, context_r)

            context_beg = max(orf_beg - 10, 0)
            context_end = orf_beg
            context = seq_local[context_beg:context_end]
            orf = seq_local[orf_beg:orf_end]
            stop_codon = seq_local[orf_end:orf_end + 3]

            orfs.append({'orf_coords': (loc[0], loc[1]),
                         'sc_score': sc_score,
                         'context': context,
                         'orf_seq': orf,
                         'stop_codon': stop_codon,
                         'frame': orf_frame})

    for orf in orfs:
        orf_start = orf['orf_coords'][0]
        orf_end = orf['orf_coords'][1]
        orf_len = orf_end - orf_start
        ovrlp = overlap((hit_start, hit_end), (orf_start, orf_end))
        ovrlp = ovrlp / max(orf_len, hit_len)
        orf['ovrlp'] = ovrlp
        # orf['grade'] = ovrlp
        orf['grade'] = (ovrlp ** 1.25) * orf['sc_score']

    orfs = sorted(orfs, key=itemgetter('grade'), reverse=True)

    good_orfs = list()
    bad_orfs = list()

    for orf in orfs:
        orf_start = orf['orf_coords'][0]
        orf_end = orf['orf_coords'][1]
        sc_score = orf['sc_score']
        context = orf['context']
        orf_seq = orf['orf_seq']
        orf_frame = orf['frame']
        ovrlp = orf['ovrlp']
        orf_grade = orf['grade']
        stop_codon = orf['stop_codon']

        tmp = orf_seq + stop_codon

        log_str += ('{:.4f}'.format(orf_grade) + ' '
                    + '{:.4f}'.format(ovrlp) + ' '
                    + '{:.4f}'.format(sc_score) + ' '
                    + str(len(tmp)).rjust(5) + ' '
                    + context.rjust(10) + ' '
                    + orf_seq[0:3] + ' '
                    + orf_seq[3:13] + ' '
                    + '...' + ' '
                    + orf_seq[-33:] + ' '
                    + stop_codon.ljust(3))

        if min_overlap > 1:
            min_overlap = min_overlap / 100

        min_overlap = max(min_overlap - 0.1, 0.75)

        orf_to_check = (orf_start, orf_end, orf_frame, orf_grade)
        passed = True

        if ovrlp < min_overlap:
            log_str += ' overlap;'
            passed = False
        if len(tmp) < min_len:
            log_str += ' min_len;'
            passed = False
        if len(tmp) > max_len:
            log_str += ' max_len;'
            passed = False
        if allow_no_strt_cod is False and tmp[0:3] not in start_codons:
            log_str += ' no_start_codon;'
            passed = False
        if allow_no_stop_cod is False and tmp[-3:] not in stop_codons:
            log_str += ' no_stop_codon;'
            passed = False
        if orf_frame != frame:
            log_str += ' frame;'
            passed = False

        if passed is True:
            log_str += ' <<<'
            good_orfs.append(orf_to_check)
        else:
            bad_orfs.append(orf_to_check)

        log_str += '\n'

    return good_orfs, bad_orfs, log_str
