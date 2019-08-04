# -*- coding: utf-8 -*-
"""Parsers."""

from kakapo.py_v_diffs import unicode
from io import StringIO
from xml.etree import ElementTree
from xmltodict import parse as parse_xml


def parse_esummary_xml_text(esummary_xml_text):  # noqa
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


def parse_efetch_sra_xml_text(efetch_sra_xml_text):  # noqa
    return [parse_xml(efetch_sra_xml_text)['SraRunInfo']['Row']]


def parse_gbseq_xml_text(gbseq_xml_text):
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
            temp_acc_ver = accession + '.' + version

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
