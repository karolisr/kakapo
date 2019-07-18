# -*- coding: utf-8 -*-
"""This module reads and writes of biological sequence and alignment files."""


from __future__ import print_function

from xml.etree import ElementTree
from datetime import datetime

from kakapo.py_v_diffs import handle_types
from kakapo.seq import SeqRecord


def _parse_gbseq_xml_handle(gbseq_xml_handle):
    """Parse GBSeq XML handle.

    :param gbseq_xml_handle: A handle to parse.
    :type gbseq_xml_handle: handle

    :returns: A dictionary with these keys: accession, date_create,
        date_update, definition, division, features, length, mol_type,
        organism, seq, strandedness, taxid, topology, version
    :rtype: dict
    """
    tree = ElementTree.parse(gbseq_xml_handle)
    root = tree.getroot()

    return_value = list()

    for rec in root.findall('GBSeq'):

        accession = None
        version = None
        # gi = None

        # n_seqids = rec.find('GBSeq_other-seqids')
        # for n_seqid in n_seqids.findall('GBSeqid'):
        #     if n_seqid.text.startswith('gi|'):
        #         gi = n_seqid.text.strip('gi|')

        temp_acc_ver = rec.find('GBSeq_accession-version')
        temp_acc_ver = temp_acc_ver.text
        accession = temp_acc_ver.split('.')[0]
        # version = int(temp_acc_ver.split('.')[1])
        version = temp_acc_ver.split('.')[1]

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

        # taxonomy = rec.find('GBSeq_taxonomy')
        # taxonomy = taxonomy.text.split('; ')

        seq = rec.find('GBSeq_sequence')
        length = rec.find('GBSeq_length')
        if seq is not None:
            seq = seq.text
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
                        start = start - 1
                    elif end < start:
                        end = end - 1
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
                        qualifiers.append(
                            {n_qualifier_name.text: n_qualifier_value.text})

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
        record_dict['mol_type'] = mol_type

        record_dict['accession'] = accession
        record_dict['version'] = version
        # record_dict['gi'] = gi
        record_dict['definition'] = definition

        record_dict['strandedness'] = strandedness
        record_dict['topology'] = topology
        record_dict['division'] = division

        record_dict['date_create'] = date_create
        record_dict['date_update'] = date_update

        record_dict['taxid'] = int(taxid)
        record_dict['organism'] = organism
        # record_dict['taxonomy'] = taxonomy

        record_dict['features'] = features

        record_dict['length'] = length

        return_value.append(record_dict)

    return return_value


def seq_records_from_gbseq_xml_handle(gbseq_xml_handle):
    parsed = _parse_gbseq_xml_handle(gbseq_xml_handle)
    return_value = list()
    for r in parsed:
        seq_record = SeqRecord(
            seq=r['seq'],
            mol_type=r['mol_type'],
            accession=r['accession'],
            version=r['version'],
            description=r['definition'],
            strandedness=r['strandedness'],
            topology=r['topology'],
            division=r['division'],
            date_create=(datetime.strptime(r['date_create'],
                                           '%d-%b-%Y')).strftime('%Y-%m-%d'),
            date_update=(datetime.strptime(r['date_update'],
                                           '%d-%b-%Y')).strftime('%Y-%m-%d'),
            taxid=r['taxid'],
            organism=r['organism'],
            features=r['features'])
        return_value.append(seq_record)
    return return_value


def seq_records_from_dict_list(seq_record_dict_list):
    return_value = list()
    for r in seq_record_dict_list:
        seq_record = SeqRecord(
            seq=r['seq'],
            mol_type=r['mol_type'],
            accession=r['accession'],
            version=r['version'],
            description=r['definition'],
            strandedness=r['strandedness'],
            topology=r['topology'],
            division=r['division'],
            date_create=(datetime.strptime(r['date_create'],
                                           '%d-%b-%Y')).strftime('%Y-%m-%d'),
            date_update=(datetime.strptime(r['date_update'],
                                           '%d-%b-%Y')).strftime('%Y-%m-%d'),
            taxid=r['taxid'],
            organism=r['organism'],
            features=r['features'])
        return_value.append(seq_record)
    return return_value


def _parse_esummary_xml_handle(esummary_xml_handle):
    tree = ElementTree.parse(esummary_xml_handle)
    root = tree.getroot()

    return_value = list()

    for rec in root.findall('DocSum'):

        record_dict = dict()

        temp_id = rec.find('Id')
        temp_id = temp_id.text

        record_dict['id'] = temp_id

        for itm in rec.findall('Item'):

            # print(itm.tag, itm.attrib, itm.text)
            # itm.attrib['Name']
            # itm.attrib['Type']

            record_dict[itm.attrib['Name']] = itm.text

        return_value.append(record_dict)

    return return_value


def write_fasta_file(records, file_path_or_handle):
    handle = False

    for h in handle_types:
        if issubclass(type(file_path_or_handle), h):
            handle = True
            break

    if not handle:
        file_path_or_handle = open(file_path_or_handle, 'w')

    for rec in records:

        description = None
        seq = None

        if hasattr(rec, 'description'):
            description = rec.description
        elif hasattr(rec, 'name'):
            description = rec.name
        elif isinstance(rec, dict):
            if 'description' in rec:
                description = rec['description']
            elif 'name' in rec:
                description = rec['name']
        else:
            raise Exception('No "name" or description attribute or key in record.')

        if hasattr(rec, 'seq'):
            seq = rec.seq.seq
        elif isinstance(rec, dict):
            if 'seq' in rec:
                seq = rec['seq']
                if hasattr(seq, 'seq'):
                    seq = seq.seq
        else:
            raise Exception('No "seq" attribute or key in record.')

        fasta_entry = '>' + description + '\n' + seq
        file_path_or_handle.write(fasta_entry + '\n')

    if not handle:
        file_path_or_handle.close()


def read_fasta_file(file_path_or_handle):
    handle = False

    for h in handle_types:
        if issubclass(type(file_path_or_handle), h):
            handle = True
            break

    if not handle:
        file_path_or_handle = open(file_path_or_handle, 'r')

    seq_names = list()
    seqs = list()
    seq_name = None
    seq = ''
    for ln in file_path_or_handle:
        ln = ln.strip('\n')
        if ln.startswith('>'):
            if seq_name is not None:
                seqs.append(seq)
                seq = ''
            seq_name = ln.strip('>')
            seq_names.append(seq_name)
        else:
            seq = seq + ln

    seqs.append(seq)

    if not handle:
        file_path_or_handle.close()

    return_list = list()
    seq_list = list(zip(seq_names, seqs))
    for s in seq_list:
        rec = {'description': s[0], 'seq': s[1].upper()}
        return_list.append(rec)

    return return_list
