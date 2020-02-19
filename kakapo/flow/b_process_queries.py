# -*- coding: utf-8 -*-

"""Kakapo workflow: Process Queries."""

import pickle
import datetime

from io import StringIO
from os import remove as osremove
from os.path import exists as ope
from os.path import join as opj

from kakapo.bioio import filter_fasta_text_by_length
from kakapo.bioio import read_fasta
from kakapo.bioio import standardize_fasta_text
from kakapo.bioio import write_fasta
from kakapo.config import PICKLE_PROTOCOL
from kakapo.ebi_domain_search import pfam_entry
from kakapo.ebi_domain_search import pfam_seqs
from kakapo.ebi_domain_search import prot_ids_for_tax_ids
from kakapo.ebi_proteins import fasta_by_accession_list
from kakapo.entrez import accessions as accessions_ncbi
from kakapo.entrez import dnld_seqs as dnld_ncbi_seqs
from kakapo.entrez import esearch
from kakapo.entrez import esummary as entrez_summary
from kakapo.seq import translate
from kakapo.seq import untranslate
from kakapo.translation_tables import TranslationTable
from kakapo.vsearch import run_cluster_fast


def pfam_uniprot_accessions(ss, pfam_acc, tax_ids, dir_cache_pfam_acc,
                            linfo=print):
    if len(pfam_acc) > 0:
        linfo('Downloading UniProt accessions for Pfam accessions [' + ss + ']')
    pfam_seqs_list = []
    for pa in pfam_acc:
        pfam_id = pfam_entry(pa)[0]['id']
        linfo(pa + ': ' + pfam_id)
        __ = opj(dir_cache_pfam_acc, pa + '__' + ss)
        if ope(__):
            with open(__, 'rb') as f:
                acc = pickle.load(f)
            pfam_seqs_list = pfam_seqs_list + acc
        else:
            # Note: the results may include "obsolete" accessions.
            # This is not a problem, they will not appear in the set of
            # downloaded sequences from UniProt.
            acc = pfam_seqs(query=pa)
            pfam_seqs_list = pfam_seqs_list + acc
            with open(__, 'wb') as f:
                pickle.dump(acc, f, protocol=PICKLE_PROTOCOL)

    pfam_uniprot_acc = prot_ids_for_tax_ids(pfam_seqs_list, tax_ids)
    return pfam_uniprot_acc


def dnld_pfam_uniprot_seqs(ss, uniprot_acc, aa_uniprot_file, dir_cache_prj,
                           linfo=print):
    if len(uniprot_acc) != 0:
        __ = opj(dir_cache_prj, 'aa_uniprot_acc_cache__' + ss)
        prev_uniprot_acc = []
        if ope(__):
            with open(__, 'rb') as f:
                prev_uniprot_acc = pickle.load(f)

        with open(__, 'wb') as f:
            pickle.dump(uniprot_acc, f, protocol=PICKLE_PROTOCOL)

        if (set(uniprot_acc) != set(prev_uniprot_acc)) or \
           (not ope(aa_uniprot_file)):

            linfo('Downloading Pfam protein sequences from UniProt [' + ss + ']')
            # Note: the number of sequences downloaded from UniProt may
            # be less than the total number of accessions. This is normal
            # as Pfam may return "obsolete" accessions, which will not be
            # downloaded here.
            __ = fasta_by_accession_list(uniprot_acc)
            __ = standardize_fasta_text(__)

            with open(aa_uniprot_file, 'w') as f:
                f.write(__)
    else:
        if ope(aa_uniprot_file):
            osremove(aa_uniprot_file)


def user_entrez_search(ss, queries, dir_cache_prj, ncbi_longevity,
                       linfo=print):
    dnld_needed = True
    accs = []
    if len(queries) != 0:

        time_stamp_now = datetime.datetime.now()
        time_stamp_file = opj(dir_cache_prj, 'ncbi_prot_time_stamp__' + ss)
        time_stamp = None
        if ope(time_stamp_file):
            with open(time_stamp_file, 'rb') as f:
                time_stamp = pickle.load(f)
                time_diff = time_stamp_now - time_stamp
                if time_diff < ncbi_longevity:
                    dnld_needed = False

        if dnld_needed is True:
            linfo('Searching for protein sequences on NCBI [' + ss + ']')
            for q in queries:
                accs = accs + accessions_ncbi(esearch(term=q, db='protein'))
            with open(time_stamp_file, 'wb') as f:
                pickle.dump(datetime.datetime.now(), f,
                            protocol=PICKLE_PROTOCOL)
        else:
            days = ncbi_longevity.total_seconds() / 60 / 60 / 24
            days = '{:.2f}'.format(days)
            linfo('NCBI results are less than ' + days + ' day(s) old. Will not search again. [' + ss + ']')
            pickle_file = opj(dir_cache_prj, 'ncbi_prot_metadata_cache__' + ss)
            if ope(pickle_file):
                with open(pickle_file, 'rb') as f:
                    pickled = pickle.load(f)
                    accs = [x['accessionversion'] for x in pickled]

    return accs


def user_protein_accessions(ss, prot_acc_user, dir_cache_prj, taxonomy,
                            linfo=print):
    if len(prot_acc_user) > 0:
        linfo('Reading user provided protein accessions [' + ss + ']')
        pickle_file = opj(dir_cache_prj, 'ncbi_prot_metadata_cache__' + ss)
        acc_old = set()
        if ope(pickle_file):
            with open(pickle_file, 'rb') as f:
                pickled = pickle.load(f)
                acc_old = set([x['accessionversion'] for x in pickled])

        if acc_old == set(prot_acc_user):
            pa_info = pickled
        else:
            pa_info = entrez_summary(prot_acc_user, 'protein')

        prot_acc = []
        prot_info_to_print = []
        for pa in pa_info:
            acc = pa['accessionversion']
            prot_acc.append(acc)
            title = pa['title']
            title_split = title.split('[')
            if len(title_split) == 2:
                organism = title_split[1].replace(']', '').strip().replace('_', ' ')
            else:
                taxid = pa['taxid']
                organism = taxonomy.scientific_name_for_taxid(taxid)
            title = title_split[0]
            title = title.lower().strip()
            title = title.replace('_', ' ').replace('-', ' ')
            title = title.replace(',', '')
            title = title[0].upper() + title[1:] + ' [' + organism + ']'
            prot_info_to_print.append((title, acc))

        prot_info_to_print = sorted(prot_info_to_print)
        for pi in prot_info_to_print:
            title = pi[0]
            acc = pi[1]
            if len(title) > 80:
                title = title[:77] + '...'
            linfo(acc + ': ' + title)

        with open(pickle_file, 'wb') as f:
            pickle.dump(pa_info, f, protocol=PICKLE_PROTOCOL)

        return prot_acc

    else:

        return prot_acc_user


def dnld_prot_seqs(ss, prot_acc_user, aa_prot_ncbi_file, linfo=print):
    if len(prot_acc_user) != 0:
        acc_old = set()
        if ope(aa_prot_ncbi_file):
            _ = read_fasta(aa_prot_ncbi_file)
            acc_old = set([x.split('|')[0] for x in tuple(_.keys())])

        if acc_old == set(prot_acc_user):
            return prot_acc_user
        else:
            linfo('Downloading protein sequences from NCBI [' + ss + ']')
            _ = dnld_ncbi_seqs(prot_acc_user, 'protein')
            prot_acc_user_new = list()
            for rec in _:
                accession = rec['accession']
                version = rec['version']
                defn = rec['definition']
                organism = rec['organism']

                if version is not None:
                    new_acc = accession + '.' + version
                else:
                    new_acc = accession

                prot_acc_user_new.append(new_acc)

                defn_new = defn.split('[' + organism + ']')[0]
                defn_new = defn_new.lower().strip()
                defn_new = defn_new.replace(' ', '_').replace('-', '_')
                defn_new = defn_new.replace(',', '')
                defn_new = defn_new[0].upper() + defn_new[1:]

                rec['definition'] = defn_new

            prot_acc_user = prot_acc_user_new
            write_fasta(_, aa_prot_ncbi_file)
    else:
        if ope(aa_prot_ncbi_file):
            osremove(aa_prot_ncbi_file)

    return prot_acc_user


def user_aa_fasta(ss, user_queries, aa_prot_user_file, linfo=print):
    _ = ''
    if len(user_queries) > 0:
        linfo('Reading user provided AA sequences [' + ss + ']')
        for ap in user_queries:
            linfo(ap)
            with open(ap, 'r') as f:
                _ = _ + f.read()
    if _ != '':
        with open(aa_prot_user_file, 'w') as f:
            f.write(standardize_fasta_text(_))


def combine_aa_fasta(ss, fasta_files, aa_queries_file, linfo=print):
    linfo('Combining all AA query sequences [' + ss + ']')
    _ = ''
    for fasta_file in fasta_files:
        if ope(fasta_file):
            with open(fasta_file, 'r') as f:
                _ = _ + f.read()

    if _ != '':
        with open(aa_queries_file, 'w') as f:
            f.write(_)
    else:
        # TODO: Do not exit here. Stop at the point in the workflow
        # when the queries are actually needed.
        linfo('No queries were provided for ' + ss + '. Exiting.')
        exit(0)


def filter_queries(ss, aa_queries_file, min_query_length,
                   max_query_length, max_query_identity, vsearch,
                   prot_acc_user, overwrite, linfo=print):
    _ = ''
    with open(aa_queries_file, 'r') as f:
        _ = f.read()

    linfo('Filtering AA query sequences [' + ss + ']')
    linfo('min_query_length: ' + str(min_query_length))
    linfo('max_query_length: ' + str(max_query_length))
    linfo('max_query_identity: ' + str(max_query_identity))

    _ = filter_fasta_text_by_length(_, min_query_length, max_query_length)

    tmp1 = aa_queries_file + '_temp1'
    tmp2 = aa_queries_file + '_temp2'
    tt = TranslationTable(1)
    parsed_fasta_1 = read_fasta(StringIO(_))
    for dfn in parsed_fasta_1:
        parsed_fasta_1[dfn] = untranslate(parsed_fasta_1[dfn], tt.table_inv)
    write_fasta(parsed_fasta_1, tmp1)
    run_cluster_fast(vsearch, max_query_identity, tmp1, tmp2)
    parsed_fasta_2 = read_fasta(tmp2)
    prot_acc_user_new = list()
    for dfn in parsed_fasta_2:
        parsed_fasta_2[dfn] = translate(parsed_fasta_2[dfn], tt.table,
                                        tt.start_codons)
        acc = dfn.split('|')[0]
        if acc in prot_acc_user:
            prot_acc_user_new.append(acc)

    if overwrite is True:
        write_fasta(parsed_fasta_2, aa_queries_file)

    osremove(tmp1)
    osremove(tmp2)

    return prot_acc_user_new
