# -*- coding: utf-8 -*-

"""
Parse kakapo project configuration (.ini) file
"""

import re
from os.path import abspath
from os.path import basename
from os.path import dirname
from os.path import expanduser
from os.path import join

from kakapo.py_v_diffs import ConfigParser
from kakapo.helpers import list_of_files
from kakapo.helpers import replace_line_in_file
from kakapo.ebi_domain_search import pfam_entry


def _parse_taxa(taxa, tax_group, taxonomy, config_file_path):
    txids = list()

    for tax in taxa:
        if tax.isdigit():
            txids.append(int(tax))
        else:
            tax_orig = tax
            txid = taxonomy.tax_id_for_name_and_group_tax_id(
                name=tax, group_tax_id=tax_group)

            if txid is None:
                msg = 'NCBI taxonomy ID for ' + tax + ' could not be found.'
                print(msg)
                replace_line_in_file(
                    file_path=config_file_path,
                    line_str=tax_orig,
                    replace_str='; NCBI taxid not found: ' + tax)

            else:
                txids.append(int(txid))
                msg = 'NCBI taxonomy ID for ' + tax + ' is ' + str(txid)
                print(msg)
                replace_line_in_file(
                    file_path=config_file_path,
                    line_str=tax_orig,
                    replace_str='; ' + tax + '\n' + str(txid))

    return txids


def _parse_pfam(pfam_entries, config_file_path):
    pfam_acc = list()

    for pf in pfam_entries:
        pf_match = re.match('PF\d+', pf, flags=re.IGNORECASE)
        if pf_match is not None:
            pfam_acc.append(pf)
        else:
            pf_orig = pf
            pf_entry = pfam_entry(pf)
            if len(pf_entry) == 0:
                replace_line_in_file(
                    file_path=config_file_path,
                    line_str=pf_orig,
                    replace_str='; Pfam accession not found: ' + pf)
            else:
                acc = pf_entry[0]['acc']
                pfam_acc.append(acc)
                replace_line_in_file(
                    file_path=config_file_path,
                    line_str=pf_orig,
                    replace_str='; ' + pf + '\n' + str(acc))

    return pfam_acc


def config_file_parse(file_path, taxonomy):  # noqa
    cfg = ConfigParser()
    cfg.optionxform = str
    cfg.read(file_path)

    # General
    project_name = cfg.get('General', 'project_name')
    email = cfg.get('General', 'email')
    output_directory = abspath(expanduser(cfg.get(
        'General', 'output_directory')))

    # Target SRA accessions
    sras = cfg.items('Target SRA accessions')
    sras = [x[0] for x in sras]

    # Target FASTQ files
    fastq_temp = cfg.items('Target FASTQ files')

    fq_pe = []
    fq_se = []

    for entry in fastq_temp:

        key = entry[0]
        val = entry[1]

        if key.startswith('pe_'):
            f_name = basename(val)
            d_path = abspath(expanduser(dirname(val)))
            pattern = re.escape(f_name).replace('\\*', '.')
            files = list_of_files(d_path)
            pe = [f for f in files if re.match(pattern, f) is not None]
            pe.sort()
            pe = [join(d_path, f) for f in pe]
            fq_pe.append(pe)

        elif key.startswith('se_'):
            se = abspath(expanduser(val))
            fq_se.append(se)

    # Target assemblies: FASTA files (DNA)
    assmbl = cfg.items('Target assemblies: FASTA files (DNA)')
    assmbl = [abspath(expanduser(x[0])) for x in assmbl]

    # Query taxonomic group
    tax_group_raw = cfg.items('Query taxonomic group')

    if len(tax_group_raw) != 1:
        raise Exception('One taxonomic group should be listed.')

    tax_group = tax_group_raw[0][0].lower()
    tax_group_name = tax_group.title()

    group_tax_ids = {'animals': 33208,
                     'archaea': 2157,
                     'bacteria': 2,
                     'fungi': 4751,
                     'plants': 33090,
                     'viruses': 10239}

    tax_group = group_tax_ids[tax_group]

    # Query taxa
    taxa_temp = cfg.items('Query taxa')
    taxa_temp = [x[0] for x in taxa_temp]
    tax_ids = _parse_taxa(taxa=taxa_temp,
                          tax_group=tax_group,
                          taxonomy=taxonomy,
                          config_file_path=file_path)

    # Query filters
    min_query_length = cfg.getint('Query filters', 'min_query_length')
    max_query_length = cfg.getint('Query filters', 'max_query_length')

    # Query Pfam families
    pfam_temp = cfg.items('Query Pfam families')
    pfam_temp = [x[0] for x in pfam_temp]
    pfam_acc = _parse_pfam(pfam_entries=pfam_temp, config_file_path=file_path)

    # Query NCBI protein and/or UniProt accessions
    prot_acc = cfg.items('Query NCBI protein and/or UniProt accessions')
    prot_acc = [x[0] for x in prot_acc]

    # Query FASTA files (Amino Acid)
    user_queries = cfg.items('Query FASTA files (Amino Acid)')
    user_queries = [abspath(expanduser(x[0])) for x in user_queries]

    # BLAST SRA/FASTQ
    blast_1_evalue = cfg.get('BLAST SRA/FASTQ', 'evalue')
    blast_1_max_target_seqs = cfg.get('BLAST SRA/FASTQ', 'max_target_seqs')
    blast_1_qcov_hsp_perc = cfg.get('BLAST SRA/FASTQ', 'qcov_hsp_perc')
    blast_1_culling_limit = cfg.get('BLAST SRA/FASTQ', 'culling_limit')

    # BLAST assemblies
    blast_2_evalue = cfg.get('BLAST assemblies', 'evalue')
    blast_2_max_target_seqs = cfg.get('BLAST assemblies', 'max_target_seqs')
    blast_2_qcov_hsp_perc = cfg.get('BLAST assemblies', 'qcov_hsp_perc')
    blast_2_culling_limit = cfg.get('BLAST assemblies', 'culling_limit')

    # ------------------------------------------------------------------------

    ret_dict = {'project_name': project_name,
                'email': email,
                'output_directory': output_directory,
                'sras': sras,
                'fq_pe': fq_pe,
                'fq_se': fq_se,
                'assmbl': assmbl,
                'min_query_length': min_query_length,
                'max_query_length': max_query_length,
                'user_queries': user_queries,
                'blast_1_culling_limit': blast_1_culling_limit,
                'blast_1_evalue': blast_1_evalue,
                'blast_1_max_target_seqs': blast_1_max_target_seqs,
                'blast_1_qcov_hsp_perc': blast_1_qcov_hsp_perc,
                'blast_2_culling_limit': blast_2_culling_limit,
                'blast_2_evalue': blast_2_evalue,
                'blast_2_max_target_seqs': blast_2_max_target_seqs,
                'blast_2_qcov_hsp_perc': blast_2_qcov_hsp_perc,
                'tax_group': tax_group,
                'tax_group_name': tax_group_name,
                'tax_ids': tax_ids,
                'pfam_acc': pfam_acc,
                'prot_acc': prot_acc}

    return ret_dict
