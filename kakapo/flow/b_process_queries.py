"""Kakapo workflow: Process Queries."""

import pickle
import datetime

from os import remove as osremove
from os.path import exists as ope
from os.path import join as opj

from kakapo.tools.bioio import filter_fasta_by_length
from kakapo.tools.bioio import read_fasta
from kakapo.tools.bioio import standardize_fasta_text
from kakapo.tools.bioio import write_fasta
from kakapo.tools.config import PICKLE_PROTOCOL
from kakapo.tools.ebi_domain_search import pfam_entry
from kakapo.tools.ebi_domain_search import pfam_seqs
from kakapo.tools.ebi_domain_search import prot_ids_for_tax_ids
from kakapo.tools.ebi_proteins import fasta_by_accession_list
from kakapo.tools.eutils import accs as accs_eutil
from kakapo.tools.eutils import search as search_eutil
from kakapo.tools.eutils import seqs as dnld_ncbi_seqs
from kakapo.tools.eutils import summary as summary_eutil
from kakapo.tools.seq import SEQ_TYPE_AA, SEQ_TYPE_DNA
from kakapo.tools.vsearch import run_cluster_fast

from kakapo.utils.logging import Log


def pfam_uniprot_accessions(ss, pfam_acc, tax_ids, dir_cache_pfam_acc):
    if len(pfam_acc) > 0:
        Log.msg_inf('Downloading UniProt accessions for Pfam accessions:', ss)
    pfam_seqs_list = []
    for pa in pfam_acc:
        pfam_id = pfam_entry(pa)[0]['id']
        Log.msg(pa + ':', pfam_id)
        _ = opj(dir_cache_pfam_acc, pa + '__' + ss)
        if ope(_):
            with open(_, 'rb') as f:
                acc = pickle.load(f)
            pfam_seqs_list = pfam_seqs_list + acc
        else:
            # Note: the results may include "obsolete" accessions.
            # This is not a problem, they will not appear in the set of
            # downloaded sequences from UniProt.
            acc = pfam_seqs(query=pa)
            pfam_seqs_list = pfam_seqs_list + acc
            with open(_, 'wb') as f:
                pickle.dump(acc, f, protocol=PICKLE_PROTOCOL)

    pfam_uniprot_acc = prot_ids_for_tax_ids(pfam_seqs_list, tax_ids)
    return pfam_uniprot_acc


def dnld_pfam_uniprot_seqs(ss, uniprot_acc, aa_uniprot_file, dir_cache_prj):
    if len(uniprot_acc) != 0:
        _ = opj(dir_cache_prj, 'aa_uniprot_acc_cache__' + ss)
        prev_uniprot_acc = []
        if ope(_):
            with open(_, 'rb') as f:
                prev_uniprot_acc = pickle.load(f)

        with open(_, 'wb') as f:
            pickle.dump(uniprot_acc, f, protocol=PICKLE_PROTOCOL)

        if (set(uniprot_acc) != set(prev_uniprot_acc)) or \
           (not ope(aa_uniprot_file)):

            Log.msg_inf('Downloading Pfam protein sequences from UniProt:', ss)
            # Note: the number of sequences downloaded from UniProt may
            # be less than the total number of accessions. This is normal
            # as Pfam may return "obsolete" accessions, which will not be
            # downloaded here.
            _ = fasta_by_accession_list(uniprot_acc)
            _ = standardize_fasta_text(_, SEQ_TYPE_AA, pfam=True)

            write_fasta(_, aa_uniprot_file)

    else:
        if ope(aa_uniprot_file):
            osremove(aa_uniprot_file)


def user_entrez_search(ss, queries, dir_cache_prj, requery_after):
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
                if time_diff < requery_after:
                    dnld_needed = False

        if dnld_needed is True:
            Log.msg_inf('Searching for protein sequences on NCBI:', ss)
            for q in queries:
                esearch_results = search_eutil(db='protein', term=q)
                accs = accs + accs_eutil(esearch_results)
            with open(time_stamp_file, 'wb') as f:
                pickle.dump(datetime.datetime.now(), f,
                            protocol=PICKLE_PROTOCOL)
        else:
            days = requery_after.total_seconds() / 60 / 60 / 24
            days = '{:.2f}'.format(days)
            Log.msg_inf('NCBI results are less than ' + days + ' day(s) old. Will not search again.', ss)
            pickle_file = opj(dir_cache_prj, 'ncbi_prot_metadata_cache__' + ss)
            if ope(pickle_file):
                with open(pickle_file, 'rb') as f:
                    pickled = pickle.load(f)
                    accs = [x['accessionversion'] for x in pickled]

    return accs


def user_protein_accessions(ss, prot_acc_user, dir_cache_prj, taxonomy):
    if len(prot_acc_user) > 0:
        Log.msg_inf('Reading user provided protein accessions:', ss)
        pickle_file = opj(dir_cache_prj, 'ncbi_prot_metadata_cache__' + ss)
        acc_old = set()
        if ope(pickle_file):
            with open(pickle_file, 'rb') as f:
                pickled = pickle.load(f)
                acc_old = set([x['accessionversion'] for x in pickled])

        if acc_old == set(prot_acc_user):
            pa_info = pickled
        else:
            pa_info = summary_eutil('protein', prot_acc_user)

        prot_acc = []
        prot_info_to_print = []
        for pa in pa_info:
            acc = pa['accessionversion']
            prot_acc.append(acc)
            title = pa['title']
            title_split = title.split('[')
            taxid = pa['taxid']
            if 'organism' in pa:
                organism = pa['organism']
            else:
                organism = taxonomy.scientific_name_for_taxid(taxid)
                pa['organism'] = organism
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
            Log.msg(acc + ':', title)

        with open(pickle_file, 'wb') as f:
            pickle.dump(pa_info, f, protocol=PICKLE_PROTOCOL)

        return prot_acc

    else:

        return prot_acc_user


def dnld_prot_seqs(ss, prot_acc_user, aa_prot_ncbi_file, dir_cache_prj):
    if len(prot_acc_user) != 0:
        acc_old = set()
        if ope(aa_prot_ncbi_file):
            _ = read_fasta(aa_prot_ncbi_file, SEQ_TYPE_AA)
            acc_old = set([x.definition.split('|')[0] for x in _])

        if acc_old == set(prot_acc_user):
            return prot_acc_user
        else:

            pickle_file = opj(dir_cache_prj, 'ncbi_prot_metadata_cache__' + ss)
            if ope(pickle_file):
                with open(pickle_file, 'rb') as f:
                    pa_info = pickle.load(f)

            print()
            Log.msg_inf('Downloading protein sequences from NCBI:', ss)
            _ = dnld_ncbi_seqs('protein', prot_acc_user, rettype='gb',
                               retmode='xml')
            prot_acc_user_new = list()
            for rec in _:
                acc_ver = rec.accession_version
                defn = rec.definition
                organism = rec.organism

                prot_acc_user_new.append(acc_ver)

                defn_new = defn.split('[' + organism + ']')[0]
                defn_new = defn_new.lower().strip()
                defn_new = defn_new.replace(' ', '_').replace('-', '_')
                defn_new = defn_new.replace(',', '')
                defn_new = defn_new[0].upper() + defn_new[1:]

                defn_new = acc_ver + '|' + defn_new + '|' + organism
                defn_new = defn_new.replace(' ', '_').replace('-', '_')

                rec.definition = defn_new

            prot_acc_user = prot_acc_user_new
            write_fasta(_, aa_prot_ncbi_file)
    else:
        if ope(aa_prot_ncbi_file):
            osremove(aa_prot_ncbi_file)

    return prot_acc_user


def user_aa_fasta(ss, user_queries, aa_prot_user_file):
    _ = ''
    if len(user_queries) > 0:
        print()
        Log.msg_inf('Reading user provided AA sequences:', ss)
        for ap in user_queries:
            Log.msg(ap)
            with open(ap, 'r') as f:
                _ = _ + f.read()
    if _ != '':
        with open(aa_prot_user_file, 'w') as f:
            write_fasta(standardize_fasta_text(_, SEQ_TYPE_AA), f)


def combine_aa_fasta(ss, fasta_files, aa_queries_file):
    Log.msg_inf('Combining all AA query sequences:', ss)
    _ = ''
    for fasta_file in fasta_files:
        if ope(fasta_file):
            with open(fasta_file, 'r') as f:
                _ = _ + f.read()

    with open(aa_queries_file, 'w') as f:
        f.write(_)


def filter_queries(ss, aa_queries_file, min_query_length,
                   max_query_length, max_query_identity, vsearch,
                   prot_acc_user, overwrite, logging=True):

    if logging is True:
        print()
        Log.msg_inf('Filtering AA query sequences:', ss)
        Log.msg('min_query_length:', str(min_query_length))
        Log.msg('max_query_length:', str(max_query_length))
        Log.msg('max_query_identity:', str(max_query_identity))

    parsed_fasta_1 = filter_fasta_by_length(aa_queries_file, SEQ_TYPE_AA,
                                            min_query_length,
                                            max_query_length)
    tmp1 = aa_queries_file + '_temp1'
    tmp2 = aa_queries_file + '_temp2'
    for rec in parsed_fasta_1:
        rec.seq.gc_code = 1
        rec.seq = rec.seq.untranslate()
    write_fasta(parsed_fasta_1, tmp1)
    run_cluster_fast(vsearch, max_query_identity, tmp1, tmp2)
    parsed_fasta_2 = read_fasta(tmp2, SEQ_TYPE_DNA, parse_def=True)
    prot_acc_user_new = list()
    for rec in parsed_fasta_2:
        rec.seq.gc_code = 1
        rec.seq = rec.seq.translate()
        acc = rec.accession_version
        if acc in prot_acc_user:
            prot_acc_user_new.append(acc)

    if overwrite is True:
        write_fasta(parsed_fasta_2, aa_queries_file, prepend_acc=True)

    osremove(tmp1)
    osremove(tmp2)

    return prot_acc_user_new
