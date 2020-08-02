"""GFF3."""

# https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
# https://useast.ensembl.org/info/website/upload/gff3.html
# http://gmod.org/wiki/GFF3

from collections import OrderedDict


def gff_template():
    t = OrderedDict()

    t['seqid'] = '.'
    t['source'] = '.'
    t['type'] = '.'
    t['start'] = '.'
    t['end'] = '.'
    t['score'] = '.'
    t['strand'] = '.'
    t['phase'] = '.'
    t['attributes'] = '.'

    return t


def gff_text(gff_dict):
    return '\t'.join(list(gff_dict.values())) + '\n'


def gff_blast_hit(json_dict):
    gff = ''

    for transcript in json_dict:
        ann = json_dict[transcript]['blast_hit']

        # BLAST Hit annotation -----------------------------------------------
        if 'blast_hit_begin' not in ann or \
           'blast_hit_end' not in ann or \
           'frame' not in ann:
            continue

        frame = ann['frame']
        entry_blast_hit = gff_template()
        blast_hit_begin = ann['blast_hit_begin']
        blast_hit_end = ann['blast_hit_end']
        entry_blast_hit['seqid'] = transcript
        entry_blast_hit['source'] = 'KAKAPO'
        entry_blast_hit['type'] = 'BLAST'
        entry_blast_hit['start'] = str(blast_hit_begin + 1)
        entry_blast_hit['end'] = str(blast_hit_end)
        entry_blast_hit['score'] = str(ann['evalue'])
        entry_blast_hit['strand'] = '+'
        entry_blast_hit['phase'] = str(0)
        query_name = ann['query_name'].replace('_', ' ')
        name = 'name=' + query_name + ';Hit frame ' + str(frame)
        entry_blast_hit['attributes'] = name + ';note=Merged tblastn hits;'
        gff = gff + gff_text(entry_blast_hit)

    return gff


def gff_orf_good(json_dict):
    gff = ''

    for transcript in json_dict:
        ann = json_dict[transcript]

        # ORF good annotation ------------------------------------------------
        if 'orfs_good' not in ann:
            continue

        orfs_good = ann['orfs_good']

        for good_orf in orfs_good:

            good = orfs_good[good_orf]

            i = good_orf.split('ORF')[1]

            orf_begin = good['orf_begin']
            orf_end = good['orf_end']
            orf_frame = good['orf_frame']
            entry_orf = gff_template()
            entry_orf['seqid'] = transcript
            entry_orf['source'] = 'KAKAPO'
            if int(i) == 1:
                entry_orf['type'] = 'CDS'
            else:
                entry_orf['type'] = 'ORF'
            entry_orf['start'] = str(orf_begin + 1)
            entry_orf['end'] = str(orf_end)
            entry_orf['score'] = str(good['orf_grade'])
            entry_orf['strand'] = '+'
            entry_orf['phase'] = str(0)
            tt = 'transl_table=' + good['orf_tt_id'] + ';transl_table_name=' + \
                 good['orf_tt_name']
            name = 'name=ORF{} frame {}'.format(i, orf_frame) + ';' + tt
            entry_orf['attributes'] = name
            gff = gff + gff_text(entry_orf)

            if 'orf_ips_ann' in good:
                for ips_ann in good['orf_ips_ann']:
                    entry_orf = gff_template()
                    entry_orf['seqid'] = transcript
                    entry_orf['source'] = 'InterProScan'
                    entry_orf['type'] = ips_ann['type']
                    entry_orf['start'] = str(ips_ann['beg'])
                    entry_orf['end'] = str(ips_ann['end'])
                    entry_orf['strand'] = '+'
                    entry_orf['attributes'] = ips_ann['attributes']
                    gff = gff + gff_text(entry_orf)

        # --------------------------------------------------------------------

    return gff


def gff_orf_bad(json_dict):
    gff = ''

    for transcript in json_dict:
        ann = json_dict[transcript]

        # ORF bad annotation -------------------------------------------------
        if 'orfs_bad' not in ann:
            continue

        orfs_bad = ann['orfs_bad']

        for bad_orf in orfs_bad:

            bad = orfs_bad[bad_orf]

            i = bad_orf.split('ORF')[1]

            orf_begin = bad['orf_begin']
            orf_end = bad['orf_end']
            orf_frame = bad['orf_frame']
            entry_orf = gff_template()
            entry_orf['seqid'] = transcript
            entry_orf['source'] = 'KAKAPO'
            entry_orf['type'] = 'ORF BAD'
            entry_orf['start'] = str(orf_begin + 1)
            entry_orf['end'] = str(orf_end)
            entry_orf['score'] = str(bad['orf_grade'])
            entry_orf['strand'] = '+'
            entry_orf['phase'] = str(0)
            tt = 'transl_table=' + bad['orf_tt_id'] + ';transl_table_name=' + \
                 bad['orf_tt_name']
            name = 'name=ORF_BAD{} frame {}'.format(i, orf_frame) + ';' + tt
            entry_orf['attributes'] = name
            gff = gff + gff_text(entry_orf)

        # --------------------------------------------------------------------

    return gff


def gff_from_json_dict(json_dict, gff_path=None):

    # gff = '##gff-version 3\n'

    gff = (gff_blast_hit(json_dict) +
           gff_orf_good(json_dict) +
           gff_orf_bad(json_dict))

    if gff_path is not None:
        with open(gff_path, 'w') as f:
            f.write(gff)

    return gff


# seqid - name of the chromosome or scaffold; chromosome names can be given
#     with or without the 'chr' prefix. Important note: the seq ID must be one
#     used within Ensembl, i.e. a standard chromosome name or an Ensembl
#     identifier such as a scaffold ID, without any additional content such as
#     species or assembly. See the example GFF output below.

# source - name of the program that generated this feature, or the data source
#     (database or project name)

# type - type of feature. Must be a term or accession from the SOFA sequence
#     ontology

# start - Start position of the feature, with sequence numbering starting at 1.

# end - End position of the feature, with sequence numbering starting at 1.

# score - A floating point value.

# strand - defined as + (forward) or - (reverse).

# phase - One of '0', '1' or '2'. '0' indicates that the first base of the
#     feature is the first base of a codon, '1' that the second base is the
#     first base of a codon, and so on..

# attributes - A semicolon-separated list of tag-value pairs, providing
#     additional information about each feature. Some of these tags are
#     predefined, e.g. ID, Name, Alias, Parent - see the GFF documentation
#     for more details.
