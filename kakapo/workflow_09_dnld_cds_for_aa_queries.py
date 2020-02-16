# -*- coding: utf-8 -*-

"""Kakapo workflow: Download CDS for NCBI protein queries."""

import pickle
import re

from os.path import exists as ope
from os.path import join as opj

from kakapo.config import PICKLE_PROTOCOL
from kakapo.entrez import cds_acc_for_prot_acc
from kakapo.entrez import dnld_cds_nt_fasta as dnld_ncbi_cds_nt_fasta
from kakapo.entrez import taxids_for_accs


def dnld_cds_for_ncbi_prot_acc(ss, prot_acc_user, prot_cds_ncbi_file, tax,
                               dir_cache_prj, linfo=print):  # noqa

    # TODO: The function downloads more CDS than are strictly required,
    #       then filters out the unneeded ones. Sometimes this causes
    #       large amounts of data to be downloaded unnecessarily.

    pickle_file = opj(dir_cache_prj, 'ncbi_prot_cds_cache__' + ss)
    acc_old = set()
    if ope(pickle_file):
        with open(pickle_file, 'rb') as f:
            pickled = pickle.load(f)
            acc_old = set(pickled[0].keys())

    if acc_old == set(prot_acc_user):
        cds_fasta = pickled[1]
        taxids = pickled[2]
    else:
        linfo('Downloading CDS for the dereplicated set of the user-provided '
              'NCBI protein accessions [' + ss + ']')
        cds_acc_dict = cds_acc_for_prot_acc(prot_acc_user)
        cds_accessions = []
        for prot_acc in cds_acc_dict:
            cds_acc = cds_acc_dict[prot_acc]
            cds_accessions.append(cds_acc)
        cds_accessions = sorted(set(cds_accessions))
        cds_fasta = dnld_ncbi_cds_nt_fasta(cds_accessions)
        taxids = taxids_for_accs(prot_acc_user, 'protein')
        with open(pickle_file, 'wb') as f:
            pickle.dump((cds_acc_dict, cds_fasta, taxids), f,
                        protocol=PICKLE_PROTOCOL)

    prot_ids_used = []
    cds_seqs_fasta_list = []

    for rec in cds_fasta:
        description = rec.split('|')[1]
        prot_id = re.findall(r'\[protein_id=(.*?)\]', rec)

        if len(prot_id) == 1:

            prot_id = prot_id[0]

            if prot_id in prot_acc_user:
                if prot_id in prot_ids_used:
                    continue

                prot_ids_used.append(prot_id)

                taxid = taxids[prot_id]
                taxon = tax.scientific_name_for_taxid(taxid)
                seq = cds_fasta[rec]
                cds_acc = re.findall(r'^(.*?)\s',
                                     description)[0].split('_cds_')[0]

                prot_name = re.findall(r'\[protein=(.*?)\]', rec)
                if len(prot_name) > 0:
                    prot_name = prot_name[0]
                else:
                    prot_name = re.findall(r'\[gene=(.*?)\]', rec)
                    if len(prot_name) > 0:
                        prot_name = prot_name[0]
                    else:
                        prot_name = ''

                prot_name = prot_name.lower().strip()
                prot_name = prot_name.replace(' ', '_').replace('-', '_')
                prot_name = prot_name.replace(',', '')
                prot_name = prot_name[0].upper() + prot_name[1:]

                defn = prot_name + '__' + prot_id + '__QUERY'
                defn = defn + ' ' + prot_name.replace('_', ' ')
                defn = defn + '; ' + taxon + ', ' + str(taxid)
                defn = defn + '; ' + prot_id + '; ' + cds_acc

                rec_new = '>' + defn + '\n' + seq
                cds_seqs_fasta_list.append(rec_new)

    cds_seqs_fasta_text = '\n'.join(cds_seqs_fasta_list)

    with open(prot_cds_ncbi_file, 'w') as f:
        f.write(cds_seqs_fasta_text)
