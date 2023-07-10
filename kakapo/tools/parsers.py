"""Parsers."""


import re
from xml.etree import ElementTree
from functools import reduce
from operator import add

from kakapo.utils.misc import list_of_files_at_path_recursive
from kakapo.tools.seq import MOL_TO_SEQ_TYPE_MAP
from kakapo.tools.seq import Seq
from kakapo.tools.seq import SeqRecord


def parse_esummary_xml_text(esummary_xml_text):
    root = ElementTree.fromstring(esummary_xml_text)

    return_value = list()

    for rec in root.findall('DocSum'):

        record_dict = dict()

        temp_id = rec.find('Id')
        temp_id = temp_id.text

        record_dict['id'] = temp_id

        for itm in rec.findall('Item'):
            record_dict[itm.attrib['Name'].lower()] = itm.text

        return_value.append(record_dict)

    return return_value


# import csv
# def parse_efetch_sra_csv_text(efetch_sra_csv_text):
#     return list(csv.DictReader(efetch_sra_csv_text.splitlines()))


def parse_efetch_sra_xml_text(efetch_sra_xml_text):
    root = ElementTree.fromstring(efetch_sra_xml_text)

    return_value = list()

    for rec in root.findall('Row'):

        record_dict = dict()

        for itm in rec:
            record_dict[itm.tag] = itm.text

        return_value.append(record_dict)

    return return_value


# ToDo: Requires careful further review.
def parse_gb_location(s):
    tmp = s

    # Assumes that all the pieces are in the same direction.
    complement = re.findall(r'complement\((.*)\)', tmp)
    if len(complement) == 1:
        complement = complement[0]
    elif len(complement) == 0:
        complement = False
    else:
        raise Exception('Found more than one \'complement\' directive.')

    if complement is not False:
        tmp = complement
        complement = True

    join = re.findall(r'join\((.*?)\)', tmp)
    if len(join) == 1:
        join = join[0]
    elif len(join) == 0:
        join = False
    else:
        raise Exception('Found more than one \'join\' directive.')

    if join is not False:
        tmp = join

    # Assumes that all the pieces refer to the same external accession.
    ext = re.findall(r'^(.*?:)', tmp)
    if len(ext) == 1:
        ext = ext[0]
        tmp = tmp.replace(ext, '')
        ext = ext.strip(':')
    elif len(ext) == 0:
        ext = None

    if re.findall(r',(.*?:)', tmp):
        raise Exception(
            'Found more than one reference to another GenBank record.')

    tmp = tmp.split(',')
    tmp = [x.split('..') for x in tmp]

    trimmed_left = False
    if tmp[0][0].startswith('<'):
        tmp[0][0] = tmp[0][0].strip('<')
        trimmed_left = True

    trimmed_right = False
    if tmp[-1][-1].startswith('>'):
        tmp[-1][-1] = tmp[-1][-1].strip('>')
        trimmed_right = True

    # print('-' * 80)
    # print(tmp)
    # print('-' * 80)

    tmp = [(int(x[0]), int(x[1])) for x in tmp]

    if complement is True:
        tmp.reverse()

    return {'location': tmp, 'rev_comp': complement,
            'trimmed_l': trimmed_left, 'trimmed_r': trimmed_right,
            'external_ref': ext}


def eutils_loc_str(loc_parsed):
    tmp = ''
    for l in loc_parsed['location']:
        tmp += str(l[0]) + ':' + str(l[1])
        if loc_parsed['rev_comp'] is True:
            tmp += ':2'
        tmp += ','

    return tmp.strip(',')


def parse_gbseq_xml_text(gbseq_xml_text: str) -> list:
    """
    Parse GBSeq XML text.

    :param gbseq_xml_text: XML text to parse.
    :type gbseq_xml_text: str

    :returns: A dictionary with these keys: accession, date_create,
        date_update, definition, division, features, length, mol_type,
        organism, seq, strandedness, taxid, topology, version
    :rtype: dict
    """

    # -----------------------------------------------
    # print(gbseq_xml_text)
    # with open('gbseq_xml_text.xml', 'w') as f:
    #     f.write(gbseq_xml_text)
    # -----------------------------------------------

    root = ElementTree.fromstring(gbseq_xml_text)

    return_value = list()

    for rec in root.findall('GBSeq'):

        version = None

        temp_acc_ver = rec.find('GBSeq_accession-version')
        temp_acc_ver = temp_acc_ver.text.split('.')

        accession = temp_acc_ver[0]

        if len(temp_acc_ver) == 2:
            version = temp_acc_ver[1]
            temp_acc_ver = accession + '.' + version

        dbsource = rec.find('GBSeq_source-db')
        if dbsource is not None:
            dbsource = dbsource.text

        strandedness = rec.find('GBSeq_strandedness')
        if strandedness is not None:
            strandedness = strandedness.text

        date_create = rec.find('GBSeq_create-date')
        date_create = date_create.text

        date_update = rec.find('GBSeq_update-date')
        date_update = date_update.text

        division = rec.find('GBSeq_division')
        division = division.text

        topology = rec.find('GBSeq_topology')
        topology = topology.text

        definition = rec.find('GBSeq_definition')
        definition = definition.text

        mol_type = rec.find('GBSeq_moltype')
        mol_type = mol_type.text

        organism = rec.find('GBSeq_organism')
        organism = organism.text

        organelle = None
        coded_by = None
        gc_id = None

        seq = rec.find('GBSeq_sequence')
        if seq is not None:
            seq = seq.text

        # hack, for now. One GB record references another: -------------------
        # KL402854 has a reference to:
        #       <GBSeq_contig>join(APNO01005115.1:1..45588)</GBSeq_contig>
        # run:
        # eutils.seqs('nuccore', ['KL402854'])
        contig = rec.find('GBSeq_contig')
        if contig is not None:
            contig = parse_gb_location(contig.text)
            external_ref_acc = contig['external_ref']
            from kakapo.tools.eutils import seqs
            seq = seqs('nucleotide', [external_ref_acc])[0].seq
        # --------------------------------------------------------------------

        length = rec.find('GBSeq_length')
        if seq is not None:
            # seq = seq.text
            length = int(length.text)
            if length != len(seq):
                message = (
                    'Reported sequence length does not match the actual '
                    'length of sequence. Reported: {r}; Actual: {a}.')
                message = message.format(r=length, a=len(seq))
                raise Exception(message)
        else:
            seq = None
            length = None

        features = list()
        n_feature_table = rec.find('GBSeq_feature-table')
        for n_feature in n_feature_table.findall('GBFeature'):
            feature = dict()
            n_feature_key = n_feature.find('GBFeature_key')
            fk = n_feature_key.text
            feature['key'] = fk
            feature['intervals'] = list()
            feature['interval_directions'] = list()
            n_intervals = n_feature.find('GBFeature_intervals')

            for n_interval in n_intervals.findall('GBInterval'):

                n_interval_from = n_interval.find('GBInterval_from')
                n_interval_to = n_interval.find('GBInterval_to')
                n_interval_point = n_interval.find('GBInterval_point')

                n_interval_accession = n_interval.find('GBInterval_accession')
                n_interval_accession = n_interval_accession.text
                # Sometimes an interval refers to an interval in another
                # accession. I choose to ignore this. Especially since, the
                # other accession is often also picked up independently.
                # e.g. X00196.1 refers to X00198.1
                if n_interval_accession != temp_acc_ver:
                    continue

                interval = None
                interval_direction = 1

                if (n_interval_from is not None) and \
                        (n_interval_to is not None):

                    start = int(n_interval_from.text)
                    end = int(n_interval_to.text)

                    # Make intervals zero-indexed.
                    if start <= end:
                        start -= 1
                    elif end < start:
                        end -= 1
                        interval_direction = -1

                    interval = [start, end]

                elif n_interval_point is not None:

                    interval = [int(n_interval_point.text) - 1]
                    interval_direction = 0

                feature['intervals'].append(interval)
                feature['interval_directions'].append(interval_direction)

            feature['qualifiers'] = list()
            qualifiers = feature['qualifiers']
            n_qualifiers = n_feature.find('GBFeature_quals')
            if n_qualifiers is not None:
                for n_qualifier in n_qualifiers.findall('GBQualifier'):
                    n_qualifier_name = n_qualifier.find('GBQualifier_name')
                    n_qualifier_value = n_qualifier.find('GBQualifier_value')
                    if n_qualifier_value is None:
                        qualifiers.append({n_qualifier_name.text: True})
                    else:
                        qualifier_value = n_qualifier_value.text
                        qualifier_name = n_qualifier_name.text

                        if fk == 'source' and qualifier_name == 'organelle':
                            organelle = qualifier_value

                        if fk == 'CDS' and qualifier_name == 'coded_by':
                            coded_by = parse_gb_location(qualifier_value)

                        if fk == 'CDS' and qualifier_name == 'transl_table':
                            gc_id = int(qualifier_value)

                        qualifiers.append({qualifier_name: qualifier_value})

            features.append(feature)

        taxid = None
        for f in features:
            if f['key'] == 'source':
                for q in f['qualifiers']:
                    if 'db_xref' in q.keys():
                        taxid_temp = q['db_xref']
                        if taxid_temp.startswith('taxon:'):
                            taxid = taxid_temp.split('taxon:')[1]
                            break

        record_dict = dict()

        record_dict['seq'] = seq
        record_dict['contig'] = contig
        record_dict['mol_type'] = mol_type

        record_dict['db_source'] = dbsource

        record_dict['accession'] = accession
        record_dict['version'] = version
        record_dict['definition'] = definition

        record_dict['strandedness'] = strandedness
        record_dict['topology'] = topology
        record_dict['division'] = division

        record_dict['date_create'] = date_create
        record_dict['date_update'] = date_update

        record_dict['taxid'] = int(taxid)
        record_dict['organism'] = organism

        record_dict['organelle'] = organelle
        record_dict['coded_by'] = coded_by
        record_dict['gc_id'] = gc_id

        record_dict['features'] = features

        record_dict['length'] = length

        return_value.append(record_dict)

    return return_value


def seq_records_gb(gbseq_xml_text: str) -> list:
    parsed = parse_gbseq_xml_text(gbseq_xml_text)

    records = list()
    for r in parsed:

        seq_str = r['seq']
        mol_type = r['mol_type']
        accession = r['accession']
        version = r['version']
        definition = r['definition']
        taxid = r['taxid']
        organism = r['organism']
        features = r['features']
        organelle = r['organelle']
        coded_by = r['coded_by']
        gc_id = r['gc_id']

        seq = Seq(seq_str, MOL_TO_SEQ_TYPE_MAP[mol_type], gc_id)
        seq_record = SeqRecord(definition, seq)
        seq_record.accession = accession
        seq_record.version = version
        seq_record.taxid = taxid
        seq_record.organism = organism
        seq_record.features = features
        seq_record.organelle = organelle
        seq_record.coded_by = coded_by

        records.append(seq_record)

    return records

