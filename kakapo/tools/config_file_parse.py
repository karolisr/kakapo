"""Parse kakapo project configuration (.ini) file."""

import datetime
import re
from collections import OrderedDict
from configparser import (ConfigParser, MissingSectionHeaderError,
                          NoOptionError, NoSectionError)
from copy import copy
from os.path import abspath, basename, dirname
from os.path import exists as ope
from os.path import expanduser, join
from sys import exit
from typing import Any, Union

from ncbi_taxonomy_local import Taxonomy

from kakapo.tools.ebi_domain_search import pfam_entry
from kakapo.tools.eutils import find_api_key
from kakapo.utils.logging import Log
from kakapo.utils.misc import list_of_files_at_path, replace_line_in_file

GROUP_TAX_IDS = {'animals': 33208,
                 'archaea': 2157,
                 'bacteria': 2,
                 'fungi': 4751,
                 'plants': 33090,
                 'viruses': 10239}


# FixMe: This needs to be done in a more elegant way.
#        Should inform the user that the closest taxonomic term should be
#        provided
def _parse_taxa(taxa, tax_group: int, taxonomy: Taxonomy, config_file_path: str):
    txids = list()

    for tax in taxa:
        if tax.isdigit():
            txids.append(int(tax))
        else:
            # tax_orig = tax
            txid = taxonomy.taxid_for_name_and_group_taxid(tax, tax_group)

            if txid is None:
                txid = taxonomy.taxid_for_name_and_group_taxid(
                    tax.split(' ')[0], tax_group)

            if txid is None:
                txid = taxonomy.taxid_for_name_and_group_taxid(
                    tax.split('_')[0], tax_group)

            if txid is None:
                txids.append(txid)
                msg = 'NCBI taxonomy ID could not be found for:'
                Log.wrn(msg, tax)
                # replace_line_in_file(
                #     file_path=config_file_path,
                #     line_str=tax_orig,
                #     replace_str='; NCBI taxid not found: ' + tax)

            else:
                txids.append(int(txid))
                msg = 'NCBI taxonomy ID for ' + tax + ' is:'
                Log.msg(msg, str(txid))
                # replace_line_in_file(
                #     file_path=config_file_path,
                #     line_str=tax_orig,
                #     replace_str='; ' + tax + '\n' + str(txid))

    return txids


def _parse_pfam(pfam_entries, config_file_path):
    pfam_acc = list()

    for pf in pfam_entries:
        pf_match = re.match('^\\s*PF\\d+$', pf, flags=re.IGNORECASE)
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


def ss_file_parse(file_path):

    cfg = ConfigParser(delimiters=('='), allow_no_value=True,
                       empty_lines_in_values=True)
    cfg.optionxform = str  # type: ignore
    cfg.SECTCRE = re.compile(r'\[\s*(?P<header>[^]]+?)\s*\]')

    try:
        cfg.read(file_path)
    except MissingSectionHeaderError:
        Log.err('Missing section header(s) in the provided "Search Strategies" file:', file_path)
        exit(1)

    required_options = {'organelle', 'min_query_length', 'max_query_length',
                        'max_query_identity', 'min_target_orf_length',
                        'max_target_orf_length'}

    ret_dict = OrderedDict()

    sections = cfg.sections()

    for s in sections:

        o = cfg.options(s)

        if not required_options <= set(o):
            missing = required_options - (required_options & set(o))
            Log.err('Missing required option(s):' + ', '.join(missing)
                    + 'for search strategy', s)
            exit(1)

        organelle = cfg[s]['organelle']

        if organelle not in ('nucleus', 'plastid', 'mitochondrion'):
            Log.err('Organelle "' + organelle + '" should be one of:',
                    'nucleus, plastid, or mitochondrion.')
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
            pfam_families = set(pfam_families.split('\n')) - {'', 'None'}
            pfam_families = _parse_pfam(pfam_entries=pfam_families,
                                        config_file_path=file_path)
            pfam_families = sorted(pfam_families)

        if cfg.has_option(s, 'ncbi_accessions_aa'):
            ncbi_accessions_aa = str(cfg[s]['ncbi_accessions_aa'])
            ncbi_accessions_aa = sorted(set(ncbi_accessions_aa.split('\n'))
                                        - {'', 'None'})

        if cfg.has_option(s, 'entrez_search_queries'):
            entrez_search_queries = str(cfg[s]['entrez_search_queries'])
            entrez_search_queries = sorted(
                set(entrez_search_queries.split('\n')) - {'', 'None'})

        if cfg.has_option(s, 'fasta_files_aa'):
            fasta_files_aa = str(cfg[s]['fasta_files_aa'])
            fasta_files_aa = set(fasta_files_aa.split('\n')) - {'', 'None'}
            fasta_files_aa = [abspath(expanduser(x)) for x in fasta_files_aa]
            fasta_files_aa = sorted(fasta_files_aa)

        section_dict = OrderedDict(
            {'organelle': organelle,
             'min_query_length': min_query_length,
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


def config_file_parse(file_path: str, taxonomy: Taxonomy,
                      err_on_missing: bool = True,
                      editor: bool = False):

    cfg = ConfigParser(delimiters=('='), allow_no_value=True)
    cfg.optionxform = str  # type: ignore

    try:
        cfg.read(file_path)
    except MissingSectionHeaderError:
        Log.err('Missing section header(s) in the provided '
                'configuration file:', file_path)
        exit(1)

    try:
        # General
        project_name = cfg.get('General', 'project_name')

        # email = cfg.get('General', 'email')
        email = 'kraman2@uic.edu'

        output_directory = cfg.get('General', 'output_directory')
        if editor is False:
            output_directory = abspath(expanduser(output_directory))

        should_use_colors = cfg.getboolean('General', 'use_colors')

        entrez_api_key = cfg.get('General', 'entrez_api_key')
        if entrez_api_key == '':
            entrez_api_key = find_api_key()
            if entrez_api_key is None:
                entrez_api_key = ''
        if entrez_api_key != '':
            print(f'\n\tUsing Entrez API Key: {entrez_api_key}\n')

        should_run_ipr = cfg.getboolean('General', 'run_inter_pro_scan')
        should_run_rcorrector = cfg.getboolean('General', 'run_rcorrector')

        # ====================================================================
        # Defaulting to Rcorrector BEFORE Trimmomatic

        # To enable user choice of RCOR/TRIM order, add this to config file:
        # rcorrector_before_trimmomatic = Yes

        # uncomment the line below:
        # rcorrector_before_trimmomatic = cfg.getboolean('General', 'rcorrector_before_trimmomatic')

        # and comment the line below:
        rcorrector_before_trimmomatic = True
        # ====================================================================

        prepend_assmbl = cfg.getboolean('General',
                                        'prepend_assembly_name_to_sequence_name')
        kraken_confidence = cfg.getfloat('General', 'kraken_2_confidence')
        requery_after = cfg.getint('General', 'requery_after')
        if editor is False:
            requery_after = datetime.timedelta(days=requery_after)

        # Target filters
        allow_non_aug = cfg.getboolean('Target filters',
                                       'allow_non_aug_start_codon')
        allow_no_strt_cod = cfg.getboolean('Target filters',
                                           'allow_missing_start_codon')
        allow_no_stop_cod = cfg.getboolean('Target filters',
                                           'allow_missing_stop_codon')

        # Query taxonomic group
        tax_group_raw = cfg.items('Query taxonomic group')

        if err_on_missing is True:
            if len(tax_group_raw) != 1:
                raise Exception('One taxonomic group should be listed.')

        tax_group = tax_group_raw[0][0].lower()
        tax_group_name = tax_group.title()

        tax_group = GROUP_TAX_IDS[tax_group]

        # Target SRA accessions
        sras = cfg.items('Target SRA accessions')
        sras = [x[0] for x in sras]

        all_tax_ids: set[int] = set()

        # Target FASTQ files
        fastq_temp = cfg.items('Target FASTQ files')

        fq_pe = []
        fq_se = []

        for entry in fastq_temp:

            key = None
            val = entry[0]

            if len(entry) == 2 and entry[1] != '':
                key = entry[0]
                val = entry[1]

            val = val.split(':')
            if len(val) == 1:
                tmp_genus_species = basename(val[0]).split('_')
                if len(tmp_genus_species) == 1:
                    genus = tmp_genus_species[0]
                    val = [genus, val[0]]
                elif len(tmp_genus_species) >= 2:
                    genus_species = tmp_genus_species[0] + ' ' + tmp_genus_species[1]
                    val = [genus_species, val[0]]
                else:
                    val = ['', val[0]]

            taxa_temp = [val[0].replace('_', ' ')]

            if taxonomy is not None:
                # FixMe: It is possible for tax_id to be None!
                #        What happens then?
                tax_id: Union[int, None] = _parse_taxa(
                    taxa=taxa_temp,
                    tax_group=tax_group,
                    taxonomy=taxonomy,
                    config_file_path=file_path)[0]
            else:
                try:
                    tax_id = int(taxa_temp[0])
                except ValueError as e:
                    Log.err(f'{taxa_temp[0]} cannot be converted to an integer.', 'Expected TaxID.')
                    exit(1)

            # See FixMe above.
            if tax_id is None:
                tax_id = tax_group

            f_name = basename(val[1])

            if key is None:
                if '*' in f_name:
                    key = 'pe_'
                else:
                    key = 'se_'

            if key.startswith('pe_'):
                # f_name = basename(val[1])

                if '*' not in f_name:
                    print()
                    Log.err('The "*" character was not found in the PE file pattern:', f_name)
                    Log.err('Stopping.', '')
                    exit(1)

                d_path = abspath(expanduser(dirname(val[1])))
                pattern = re.escape(f_name).replace('\\*', '.') + '$'
                files, err = list_of_files_at_path(d_path)

                if files is None:
                    assert err is not None
                    Log.err('Stopping.', err)
                    exit(1)

                if editor is False:
                    assert files is not None
                    pe = [f for f in files if re.match(pattern, basename(f)) is not None]
                    pe.sort()

                    if len(pe) == 0:
                        print()
                        Log.err('No FASTQ files match the pattern:', f_name)
                        Log.err('Stopping.', '')
                        exit(1)

                    elif len(pe) == 1:
                        print()
                        Log.err('Only one FASTQ file matches the PE file pattern:', f_name)
                        Log.err('Stopping.', '')
                        exit(1)

                    elif len(pe) > 2:
                        print()
                        Log.err(str(len(pe)) + ' FASTQ files match the PE file pattern:', f_name)
                        Log.err('Stopping.', '')
                        exit(1)

                    pe = [join(d_path, f) for f in pe]

                    for _ in pe:
                        if not ope(_):
                            print()
                            Log.err('Cannot find the FASTQ file:', _)
                            Log.err('Stopping.', '')
                            exit(1)

                    fq_pe.append([tax_id, pe])
                else:
                    pe = val[1]
                    fq_pe.append([tax_id, pe])

            elif key.startswith('se_'):
                # f_name = basename(val[1])
                d_path = abspath(expanduser(dirname(val[1])))
                files, err = list_of_files_at_path(d_path)

                if err is not None:
                    Log.err(err, '')
                    if err_on_missing is True:
                        Log.err('Stopping.', '')
                        exit(1)

                if editor is False:
                    se = abspath(expanduser(val[1]))
                    if not ope(se):
                        Log.err('Cannot find the FASTQ file:', se)
                        Log.err('Stopping.', '')
                        exit(1)
                    fq_se.append([tax_id, se])

                else:
                    se = val[1]
                    fq_se.append([tax_id, se])

            # See FixMe above.
            if tax_id != tax_group:
                all_tax_ids.add(tax_id)

        # Target assemblies: FASTA files (DNA)
        _ = cfg.items('Target assemblies: FASTA files (DNA)')
        _ = [x[0].split(':') for x in _]

        for i, val in enumerate(copy(_)):
            if len(val) == 1:
                tmp_genus_species = basename(val[0]).split('_')
                if len(tmp_genus_species) == 1:
                    genus = tmp_genus_species[0]
                    _[i] = [genus, val[0]]
                elif len(tmp_genus_species) >= 2:
                    genus_species = tmp_genus_species[0] + ' ' + tmp_genus_species[1]
                    _[i] = [genus_species, val[0]]
                else:
                    _[i] = ['', val[0]]

        taxa_temp = [x[0] for x in _]
        taxa_temp = [x.split('.')[0] for x in taxa_temp]

        if taxonomy is not None:
            # FixMe: It is possible for one of the tax_ids to be None!
            #        What happens then?
            tax_ids_temp: list[Union[int, None]] = _parse_taxa(
                taxa=taxa_temp,
                tax_group=tax_group,
                taxonomy=taxonomy,
                config_file_path=file_path)
        else:
            try:
                tax_ids_temp = [int(x) for x in taxa_temp]
            except ValueError as e:
                Log.err(f'Values: {taxa_temp} cannot be converted to an integers.', 'Expected TaxIDs.')
                exit(1)

        if editor is False:
            _ = [abspath(expanduser(x[1])) for x in _]
        else:
            _ = [x[1] for x in _]

        _ = list(zip(tax_ids_temp, _))
        assmbl_temp = [[x[0], x[1]] for x in _]

        assmbl = list()
        tax_ids: list[int] = list()
        all_assemblies_found = True
        for i, a in enumerate(assmbl_temp):
            # See FixMe above.
            if a[0] is None:
                a[0] = tax_group
            tax_ids.append(a[0])
            a_path = a[1]
            if not ope(a_path):
                Log.err('Cannot find the assembly file:', a_path)
                if editor is False:
                    all_assemblies_found = False
            assmbl.append(tuple(a))

        if all_assemblies_found is False:
            Log.err('Stopping.', '')
            exit(1)

        for tax_id in tax_ids:
            # See FixMe above.
            if tax_id != tax_group:
                all_tax_ids.add(tax_id)

        # Bowtie2 filter order
        bt2_sctn = 'Bowtie2 filter order'
        bt2_order = []
        if cfg.has_section(bt2_sctn):
            bt2_order = OrderedDict(cfg.items(bt2_sctn))

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
        Log.err('Missing required section "' + err.section
                + '" in configuration file:', file_path)
        exit(1)

    except NoOptionError as err:
        Log.err('Missing required option "' + err.option
                + '" under section "' + err.section + '" in configuration file:',
                file_path)
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
                'requery_after': requery_after,
                'entrez_api_key': entrez_api_key,
                'fq_pe': fq_pe,
                'fq_se': fq_se,
                'should_use_colors': should_use_colors,
                'should_run_rcorrector': should_run_rcorrector,
                'rcorrector_before_trimmomatic': rcorrector_before_trimmomatic,
                'should_run_ipr': should_run_ipr,
                'bt2_order': bt2_order,
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


def use_colors(file_path):

    cfg = ConfigParser(delimiters=('=',), allow_no_value=True)
    cfg.optionxform = str  # type: ignore

    try:
        cfg.read(file_path)
    except MissingSectionHeaderError:
        return False

    try:
        return cfg.getboolean('General', 'use_colors')

    except NoSectionError as err:
        return False

    except NoOptionError as err:
        return False
