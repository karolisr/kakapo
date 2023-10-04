"""Kakapo workflow: Download CDS for NCBI protein queries."""

import pickle
from os.path import exists as ope
from os.path import join as opj

from kakapo.tools.bioio import seq_records_to_dict, write_fasta
from kakapo.tools.config import PICKLE_PROTOCOL
from kakapo.tools.eutils import cds as cds_for_prot
from kakapo.utils.logging import Log


def dnld_cds_for_ncbi_prot_acc(ss, prot_acc_user, prot_cds_ncbi_file, tax,
                               dir_cache_prj):

    pickle_file = opj(dir_cache_prj, 'ncbi_prot_cds_cache__' + ss)
    acc_old = set()
    pickled = None
    if ope(pickle_file):
        with open(pickle_file, 'rb') as f:
            pickled = pickle.load(f)
            acc_old = set(pickled[0])

    if acc_old == set(prot_acc_user):
        assert pickled is not None
        cds_rec_dict = pickled[1]
        Log.log(inf='The CDS for the dereplicated set of the user-provided '
                    'NCBI protein accessions have already been '
                    'downloaded:', s=ss)
    else:
        Log.log(inf='Downloading CDS for the dereplicated set of the user-provided '
                    'NCBI protein accessions:', s=ss)
        cds_rec_dict = seq_records_to_dict(cds_for_prot(prot_acc_user),
                                           prepend_acc=True)
        with open(pickle_file, 'wb') as f:
            pickle.dump((prot_acc_user, cds_rec_dict), f,
                        protocol=PICKLE_PROTOCOL)

    write_fasta(cds_rec_dict, prot_cds_ncbi_file)
