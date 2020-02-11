# -*- coding: utf-8 -*-

"""
Parse kakapo project configuration (.ini) file
"""

import re

from collections import OrderedDict
from configparser import ConfigParser
from copy import copy
from os.path import abspath
from os.path import basename
from os.path import dirname
from os.path import exists as ope
from os.path import expanduser
from os.path import join
from sys import exit

from configparser import MissingSectionHeaderError
from configparser import NoSectionError
from configparser import NoOptionError

from kakapo.ebi_domain_search import pfam_entry
from kakapo.helpers import list_of_files
from kakapo.helpers import replace_line_in_file


def _parse_taxa(taxa, tax_group, taxonomy, config_file_path, linfo=print):
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
                linfo(msg)
                # replace_line_in_file(
                #     file_path=config_file_path,
                #     line_str=tax_orig,
                #     replace_str='; NCBI taxid not found: ' + tax)

            else:
                txids.append(int(txid))
                msg = 'NCBI taxonomy ID for ' + tax + ' is ' + str(txid)
                linfo(msg)
                # replace_line_in_file(
                #     file_path=config_file_path,
                #     line_str=tax_orig,
                #     replace_str='; ' + tax + '\n' + str(txid))

    return txids


def _parse_pfam(pfam_entries, config_file_path):
    pfam_acc = list()

    for pf in pfam_entries:
        pf_match = re.match('^\s*PF\d+$', pf, flags=re.IGNORECASE)
        if pf_match is not None:
            pfam_acc.append(pf)
        else:
            pf_orig = pf
            pf_entry = pfam_entry(pf)
            if len(pf_entry) == 0:
                replace_line_in_file(
                    file_path=config_file_path,
                    line_str=pf_orig,
                    replace_str='    ; Pfam accession not found: ' + pf)
            else:
                acc = pf_entry[0]['acc']
                pfam_acc.append(acc)
                replace_line_in_file(
                    file_path=config_file_path,
                    line_str=pf_orig,
                    replace_str='    ; ' + pf + '\n    ' + str(acc))

    return pfam_acc


def ss_file_parse(file_path, linfo=print):  # noqa

    linfo('Reading search strategies file: ' + file_path)

    cfg = ConfigParser(delimiters=('='), allow_no_value=True,
                       empty_lines_in_values=True)
    cfg.optionxform = str
    cfg.SECTCRE = re.compile(r'\[\s*(?P<header>[^]]+?)\s*\]')

    try:
        cfg.read(file_path)
    except MissingSectionHeaderError:
        linfo('Error: Missing section header(s) in the provided "Search Strategies" file: ' +
              file_path)
        exit(1)

    required_options = set(('min_query_length',
                            'max_query_length',
                            'max_query_identity',
                            'min_target_orf_length',
                            'max_target_orf_length'))

    ret_dict = OrderedDict()

    sections = cfg.sections()

    for s in sections:

        o = cfg.options(s)

        if not required_options <= set(o):
            missing = required_options - (required_options & set(o))
            print('Missing required option(s):', ', '.join(missing),
                  'for search strategy', s + '.')
            exit(1)

        min_query_length = int(cfg[s]['min_query_length'])
        max_query_length = int(cfg[s]['max_query_length'])
        max_query_identity = float(cfg[s]['max_query_identity'])
        min_target_orf_length = int(cfg[s]['min_target_orf_length'])
        max_target_orf_length = int(cfg[s]['max_target_orf_length'])

        evalue = None
        max_hsps = None
        qcov_hsp_perc = None
        best_hit_overhang = None
        best_hit_score_edge = None
        max_target_seqs = None

        pfam_families = None
        ncbi_accessions_aa = None
        entrez_search_queries = None
        fasta_files_aa = None

        if cfg.has_option(s, 'evalue'):
            evalue = float(cfg[s]['evalue'])

        if cfg.has_option(s, 'max_hsps'):
            max_hsps = int(cfg[s]['max_hsps'])

        if cfg.has_option(s, 'qcov_hsp_perc'):
            qcov_hsp_perc = float(cfg[s]['qcov_hsp_perc'])

        if cfg.has_option(s, 'best_hit_overhang'):
            best_hit_overhang = float(cfg[s]['best_hit_overhang'])

        if cfg.has_option(s, 'best_hit_score_edge'):
            best_hit_score_edge = float(cfg[s]['best_hit_score_edge'])

        if cfg.has_option(s, 'max_target_seqs'):
            max_target_seqs = int(cfg[s]['max_target_seqs'])

        if cfg.has_option(s, 'pfam_families'):
            pfam_families = str(cfg[s]['pfam_families'])
            pfam_families = set(pfam_families.split('\n')) - \
                set(('', 'None'))
            pfam_families = _parse_pfam(pfam_entries=pfam_families,
                                        config_file_path=file_path)
            pfam_families = sorted(pfam_families)

        if cfg.has_option(s, 'ncbi_accessions_aa'):
            ncbi_accessions_aa = str(cfg[s]['ncbi_accessions_aa'])
            ncbi_accessions_aa = sorted(set(ncbi_accessions_aa.split('\n')) -
                                        set(('', 'None')))

        if cfg.has_option(s, 'entrez_search_queries'):
            entrez_search_queries = str(cfg[s]['entrez_search_queries'])
            entrez_search_queries = sorted(
                set(entrez_search_queries.split('\n')) - set(('', 'None')))

        if cfg.has_option(s, 'fasta_files_aa'):
            fasta_files_aa = str(cfg[s]['fasta_files_aa'])
            fasta_files_aa = set(fasta_files_aa.split('\n')) - \
                set(('', 'None'))
            fasta_files_aa = [abspath(expanduser(x)) for x in fasta_files_aa]
            fasta_files_aa = sorted(fasta_files_aa)

        section_dict = OrderedDict(
            {'min_query_length': min_query_length,
             'max_query_length': max_query_length,
             'max_query_identity': max_query_identity,
             'min_target_orf_length': min_target_orf_length,
             'max_target_orf_length': max_target_orf_length,
             'blast_2_evalue': evalue,
             'blast_2_max_hsps': max_hsps,
             'blast_2_qcov_hsp_perc': qcov_hsp_perc,
             'blast_2_best_hit_overhang': best_hit_overhang,
             'blast_2_best_hit_score_edge': best_hit_score_edge,
             'blast_2_max_target_seqs': max_target_seqs,
             'pfam_families': pfam_families,
             'ncbi_accessions_aa': ncbi_accessions_aa,
             'entrez_search_queries': entrez_search_queries,
             'fasta_files_aa': fasta_files_aa})

        ret_dict[s] = section_dict

    return ret_dict


def config_file_parse(file_path, taxonomy, linfo=print):  # noqa

    linfo('Reading configuration file: ' + file_path)

    cfg = ConfigParser(delimiters=('='), allow_no_value=True)
    cfg.optionxform = str

    try:
        cfg.read(file_path)
    except MissingSectionHeaderError:
        linfo('Error: Missing section header(s) in the provided configuration file: ' + file_path)
        exit(1)

    try:
        # General
        project_name = cfg.get('General', 'project_name')
        email = cfg.get('General', 'email')
        output_directory = abspath(expanduser(cfg.get(
            'General', 'output_directory')))
        inter_pro_scan = cfg.getboolean('General', 'run_inter_pro_scan')
        prepend_assmbl = cfg.getboolean('General',
                                        'prepend_assembly_name_to_sequence_name')
        kraken_confidence = cfg.getfloat('General', 'kraken_2_confidence')

        # Target filters
        allow_non_aug = cfg.getboolean('Target filters',
                                       'allow_non_aug_start_codon')
        allow_no_strt_cod = cfg.getboolean('Target filters',
                                           'allow_missing_start_codon')
        allow_no_stop_cod = cfg.getboolean('Target filters',
                                           'allow_missing_stop_codon')

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

        # Target SRA accessions
        sras = cfg.items('Target SRA accessions')
        sras = [x[0] for x in sras]

        all_tax_ids = set()

        # Target FASTQ files
        fastq_temp = cfg.items('Target FASTQ files')

        fq_pe = []
        fq_se = []

        for entry in fastq_temp:

            key = entry[0]
            val = entry[1]

            val = val.split(':')
            if len(val) == 1:
                genus = basename(val[0]).split('_')[0]
                val = [genus, val[0]]

            taxa_temp = [val[0]]

            tax_ids = _parse_taxa(taxa=taxa_temp,
                                  tax_group=tax_group,
                                  taxonomy=taxonomy,
                                  config_file_path=file_path,
                                  linfo=linfo)

            if key.startswith('pe_'):
                f_name = basename(val[1])
                d_path = abspath(expanduser(dirname(val[1])))
                pattern = re.escape(f_name).replace('\\*', '.')
                try:
                    files = list_of_files(d_path, linfo=linfo)
                except Exception:
                    exit(1)
                pe = [f for f in files if re.match(pattern, f) is not None]
                pe.sort()
                pe = [join(d_path, f) for f in pe]
                fq_pe.append([tax_ids[0], pe])

            elif key.startswith('se_'):
                se = abspath(expanduser(val[1]))
                fq_se.append([tax_ids[0], se])

            all_tax_ids.add(tax_ids[0])

        # Target assemblies: FASTA files (DNA)
        assmbl_temp = cfg.items('Target assemblies: FASTA files (DNA)')
        assmbl_temp = [x[0].split(':') for x in assmbl_temp]

        for i, val in enumerate(copy(assmbl_temp)):
            if len(val) == 1:
                genus = basename(val[0]).split('_')[0]
                assmbl_temp[i] = [genus, val[0]]

        taxa_temp = [x[0] for x in assmbl_temp]
        taxa_temp = [x.split('.')[0] for x in taxa_temp]

        tax_ids = _parse_taxa(taxa=taxa_temp,
                              tax_group=tax_group,
                              taxonomy=taxonomy,
                              config_file_path=file_path,
                              linfo=linfo)

        assmbl = [abspath(expanduser(x[1])) for x in assmbl_temp]
        assmbl = list(zip(tax_ids, assmbl))

        all_assemblies_found = True
        for a in assmbl:
            a_path = a[1]
            if not ope(a_path):
                linfo('Cannot find the assembly file: ' + a_path)
                all_assemblies_found = False

        if all_assemblies_found is False:
            linfo('Stopping.')
            exit(1)

        for tax_id in tax_ids:
            all_tax_ids.add(tax_id)
        all_tax_ids = tuple(sorted(all_tax_ids))

        # Kraken2 filter order
        krkn_sctn = 'Kraken2 filter order'
        krkn_order = []
        if cfg.has_section(krkn_sctn):
            krkn_order = cfg.items(krkn_sctn)

        # BLAST SRA/FASTQ
        blast_1_evalue = cfg.getfloat('BLAST SRA/FASTQ', 'evalue')
        blast_1_max_hsps = cfg.getint('BLAST SRA/FASTQ', 'max_hsps')
        blast_1_qcov_hsp_perc = cfg.getfloat('BLAST SRA/FASTQ', 'qcov_hsp_perc')
        blast_1_best_hit_overhang = cfg.getfloat('BLAST SRA/FASTQ', 'best_hit_overhang')
        blast_1_best_hit_score_edge = cfg.getfloat('BLAST SRA/FASTQ', 'best_hit_score_edge')
        blast_1_max_target_seqs = cfg.getint('BLAST SRA/FASTQ', 'max_target_seqs')

        # BLAST assemblies
        blast_2_evalue = cfg.getfloat('BLAST assemblies', 'evalue')
        blast_2_max_hsps = cfg.getint('BLAST assemblies', 'max_hsps')
        blast_2_qcov_hsp_perc = cfg.getfloat('BLAST assemblies', 'qcov_hsp_perc')
        blast_2_best_hit_overhang = cfg.getfloat('BLAST assemblies', 'best_hit_overhang')
        blast_2_best_hit_score_edge = cfg.getfloat('BLAST assemblies', 'best_hit_score_edge')
        blast_2_max_target_seqs = cfg.getint('BLAST assemblies', 'max_target_seqs')

    except NoSectionError as err:
        linfo('Error: Missing required section "' + err.section + '" in configuration file: ' + file_path)
        exit(1)

    except NoOptionError as err:
        linfo('Error: Missing required option "' + err.option +
              '" under section "' + err.section + '" in configuration file: ' + file_path)
        exit(1)

    # ------------------------------------------------------------------------

    ret_dict = {'allow_no_stop_cod': allow_no_stop_cod,
                'allow_no_strt_cod': allow_no_strt_cod,
                'allow_non_aug': allow_non_aug,
                'assmbl': assmbl,

                'blast_1_evalue': blast_1_evalue,
                'blast_1_max_hsps': blast_1_max_hsps,
                'blast_1_qcov_hsp_perc': blast_1_qcov_hsp_perc,
                'blast_1_best_hit_overhang': blast_1_best_hit_overhang,
                'blast_1_best_hit_score_edge': blast_1_best_hit_score_edge,
                'blast_1_max_target_seqs': blast_1_max_target_seqs,

                'blast_2_evalue': blast_2_evalue,
                'blast_2_max_hsps': blast_2_max_hsps,
                'blast_2_qcov_hsp_perc': blast_2_qcov_hsp_perc,
                'blast_2_best_hit_overhang': blast_2_best_hit_overhang,
                'blast_2_best_hit_score_edge': blast_2_best_hit_score_edge,
                'blast_2_max_target_seqs': blast_2_max_target_seqs,

                'email': email,
                'fq_pe': fq_pe,
                'fq_se': fq_se,
                'inter_pro_scan': inter_pro_scan,
                'kraken_confidence': kraken_confidence,
                'krkn_order': krkn_order,
                'output_directory': output_directory,
                'prepend_assmbl': prepend_assmbl,
                'project_name': project_name,
                'sras': sras,
                'tax_group': tax_group,
                'tax_group_name': tax_group_name,
                'tax_ids': all_tax_ids
                }

    return ret_dict
