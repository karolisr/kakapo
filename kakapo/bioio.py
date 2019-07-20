# -*- coding: utf-8 -*-
"""This module reads and writes of biological sequence and alignment files."""


from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import re
from collections import OrderedDict
from datetime import datetime
from io import StringIO
from xml.etree import ElementTree

from kakapo.py_v_diffs import handle_types
from kakapo.py_v_diffs import unicode
from kakapo.seq import SeqRecord
from kakapo.entrez import esearch, epost, efetch, esummary
from xmltodict import parse as parse_xml


def _parse_gbseq_xml_text(gbseq_xml_text):
    """

    Parse GBSeq XML text.

    :param gbseq_xml_text: XML text to parse.
    :type gbseq_xml_text: str

    :returns: A dictionary with these keys: accession, date_create,
        date_update, definition, division, features, length, mol_type,
        organism, seq, strandedness, taxid, topology, version
    :rtype: dict
    """
    tree = ElementTree.parse(StringIO(unicode(gbseq_xml_text)))
    root = tree.getroot()

    return_value = list()

    for rec in root.findall('GBSeq'):

        accession = None
        version = None

        temp_acc_ver = rec.find('GBSeq_accession-version')
        temp_acc_ver = temp_acc_ver.text.split('.')

        accession = temp_acc_ver[0]

        if len(temp_acc_ver) == 2:
            version = temp_acc_ver[1]

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
        record_dict['definition'] = definition

        record_dict['strandedness'] = strandedness
        record_dict['topology'] = topology
        record_dict['division'] = division

        record_dict['date_create'] = date_create
        record_dict['date_update'] = date_update

        record_dict['taxid'] = int(taxid)
        record_dict['organism'] = organism

        record_dict['features'] = features

        record_dict['length'] = length

        return_value.append(record_dict)

    return return_value


def seq_records_from_efetch_results(efetch_results):  # noqa

    return_value = []

    for r in efetch_results:

        seq_record = SeqRecord(
            seq=r['seq'],
            mol_type=r['mol_type'],
            accession=r['accession'],
            version=r['version'],
            description=r['definition'],
            strandedness=r['strandedness'],
            topology=r['topology'],
            division=r['division'],
            date_create=(datetime.strptime(
                r['date_create'], '%d-%b-%Y')).strftime('%Y-%m-%d'),
            date_update=(datetime.strptime(
                r['date_update'], '%d-%b-%Y')).strftime('%Y-%m-%d'),
            taxid=r['taxid'],
            organism=r['organism'],
            features=r['features']
        )

        return_value.append(seq_record)

    return return_value


def _parse_esummary_xml_text(esummary_xml_text):
    tree = ElementTree.parse(StringIO(unicode(esummary_xml_text)))
    root = tree.getroot()

    return_value = list()

    for rec in root.findall('DocSum'):

        record_dict = dict()

        temp_id = rec.find('Id')
        temp_id = temp_id.text

        record_dict['id'] = temp_id

        for itm in rec.findall('Item'):
            record_dict[itm.attrib['Name']] = itm.text

        return_value.append(record_dict)

    return return_value


def esearch_epost(term, db):  # noqa
    if type(term) in [list, tuple]:
        term = ' OR '.join(term)
    esearch_results = esearch(db, term)
    id_list = esearch_results['IdList']
    epost_results = epost(db, id_list)
    return epost_results


def dnld_ncbi_seqs(term, db):  # noqa
    epost_results = esearch_epost(term, db)
    efetch_results = efetch(epost_results, _parse_gbseq_xml_text, 'gb')
    seq_records = seq_records_from_efetch_results(efetch_results)
    return seq_records


def entrez_summary(term, db):  # noqa
    epost_results = esearch_epost(term, db)
    esummary_results = esummary(epost_results, _parse_esummary_xml_text)
    return esummary_results


def _parse_efetch_sra_xml_text(efetch_sra_xml_text):
    return [parse_xml(efetch_sra_xml_text)['SraRunInfo']['Row']]


def sra_info(term):  # noqa
    epost_results = esearch_epost(term, 'sra')
    efetch_results = efetch(epost_results, _parse_efetch_sra_xml_text,
                            'runinfo')
    return efetch_results


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
            elif 'definition' in rec:
                description = rec['definition']
            elif 'name' in rec:
                description = rec['name']
        else:
            raise Exception('No "name" or description attribute or key in'
                            'record.')

        if hasattr(rec, 'accession'):
            description = rec.accession + '|' + description
        elif isinstance(rec, dict):
            if 'accession' in rec:
                description = rec['accession'] + '|' + description

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


def read_fasta_file_dict(file_path_or_handle):  # noqa
    records = read_fasta_file(file_path_or_handle)
    records = {r['description']: r['seq'] for r in records}
    return records


def parse_fasta_text(text): # noqa
    desc_lines = re.findall('\>.*', text)
    lines = text.split('\n')
    data = OrderedDict()
    for l in lines:
        if l in desc_lines:
            desc = l.strip('>')
            data[desc] = ''
        else:
            data[desc] = data[desc] + l.upper()

    return data


def standardize_fasta_text(text): # noqa
    parsed_fasta = parse_fasta_text(text)
    t = ''
    for k in parsed_fasta:
        desc = k.replace(' ', '_')
        desc = '>' + desc
        seq = parsed_fasta[k]
        t = t + desc + '\n' + seq + '\n'
    return t


def filter_fasta_text_by_length(fasta_text, min_len, max_len): # noqa
    parsed_fasta = parse_fasta_text(fasta_text)

    filtered_text = ''

    for k in parsed_fasta:
        desc = '>' + k
        seq = parsed_fasta[k]
        if len(seq) < min_len:
            continue
        if len(seq) > max_len:
            continue

        filtered_text = filtered_text + desc + '\n' + seq + '\n'

    return filtered_text
