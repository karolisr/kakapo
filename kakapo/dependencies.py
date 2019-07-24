# -*- coding: utf-8 -*-

"""
Dependencies.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import with_statement

import os
import sys
import tarfile
import zipfile

from kakapo.config import DIR_DEP, DIR_CFG, OS_ID, DIST_ID
from kakapo.helpers import list_of_dirs
from kakapo.http import download_file
from kakapo.os_diffs import DEBIAN_DISTS, REDHAT_DISTS
from kakapo.shell import call


def sra_toolkit_dir_name(path): # noqa
    ld = list_of_dirs(path=path)
    pattern = 'sratoolkit'
    dl = [d for d in ld if pattern in d]
    if len(dl) > 0:
        return dl[0]
    else:
        return ''


def blast_dir_name(path): # noqa
    ld = list_of_dirs(path=path)
    pattern = 'ncbi-blast'
    dl = [d for d in ld if pattern in d]
    if len(dl) > 0:
        return dl[0]
    else:
        return ''


def spades_dir_name(path): # noqa
    ld = list_of_dirs(path=path)
    pattern = 'SPAdes'
    dl = [d for d in ld if pattern in d]
    if len(dl) > 0:
        return dl[0]
    else:
        return ''


# Seqtk
def dep_check_seqtk(): # noqa
    url = 'https://github.com/lh3/seqtk/archive/master.zip'
    dnld_path = os.path.join(DIR_DEP, 'seqtk.zip')
    dir_bin = os.path.join(DIR_DEP, 'seqtk-master')

    try:
        seqtk = 'seqtk'
        call(seqtk)
    except Exception:
        try:
            seqtk = os.path.join(dir_bin, 'seqtk')
            call(seqtk)
        except Exception:
            print('\tSeqtk was not found on this system, trying to download.')
            download_file(url, dnld_path)
            zip_ref = zipfile.ZipFile(dnld_path, 'r')
            zip_ref.extractall(DIR_DEP)
            zip_ref.close()
            try:
                print('\tCompiling Seqtk.')
                call('make', cwd=dir_bin)
                seqtk = os.path.join(dir_bin, 'seqtk')
                call(seqtk)
            except Exception:
                print('\t\tSomething went wrong while trying to compile '
                      'Seqtk.')
                print('\t\tTry downloading and installing it manually from:')
                print('\t\thttps://github.com/lh3/seqtk')
                sys.exit(1)

    print('\tSeqtk is available:')
    print('\t\t' + seqtk + '\n')

    return seqtk


def _write_trimmomatic_adapters_file():

    path_adapters = os.path.join(DIR_CFG, 'trimmomatic_adapters.fasta')

    adapters = ('>TruSeq2_SE'
                'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG'
                '>TruSeq2_PE_f'
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
                '>TruSeq2_PE_r'
                'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG'
                '>TruSeq3_IndexedAdapter'
                'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
                '>TruSeq3_UniversalAdapter'
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
                '>PrefixPE/1'
                'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
                '>PrefixPE/2'
                'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
                '>PCR_Primer1'
                'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
                '>PCR_Primer1_rc'
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
                '>PCR_Primer2'
                'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT'
                '>PCR_Primer2_rc'
                'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG'
                '>FlowCell1'
                'TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC'
                '>FlowCell2'
                'TTTTTTTTTTCAAGCAGAAGACGGCATACGA'
                '>PrefixPE/1'
                'TACACTCTTTCCCTACACGACGCTCTTCCGATCT'
                '>PrefixPE/2'
                'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
                '>PE1'
                'TACACTCTTTCCCTACACGACGCTCTTCCGATCT'
                '>PE1_rc'
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
                '>PE2'
                'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
                '>PE2_rc'
                'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC')

    if not os.path.exists(path_adapters):
        print('\tWriting Trimmomatic adapter files:\n\t\t' +
              path_adapters + '\n')
        with open(path_adapters, mode='w') as f:
            f.write(adapters)

    return path_adapters

# Trimmomatic
def dep_check_trimmomatic(): # noqa
    url = ('http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/'
           'Trimmomatic-0.39.zip')
    dnld_path = os.path.join(DIR_DEP, 'Trimmomatic-0.39.zip')
    dir_bin = os.path.join(DIR_DEP, 'Trimmomatic-0.39')
    trimmomatic = os.path.join(dir_bin, 'trimmomatic-0.39.jar')

    if not os.path.exists(trimmomatic):
        download_file(url, dnld_path)
        zip_ref = zipfile.ZipFile(dnld_path, 'r')
        zip_ref.extractall(DIR_DEP)
        zip_ref.close()

    if not os.path.exists(trimmomatic):
        print('\tCould not download Trimmomatic.')
        sys.exit(1)

    print('\tTrimmomatic is available:')
    print('\t\t' + trimmomatic + '\n')
    path_adapters = _write_trimmomatic_adapters_file()

    return trimmomatic, path_adapters


# SRA Toolkit
def dep_check_sra_toolkit(): # noqa
    if OS_ID == 'mac':
        url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/'
               'sratoolkit.current-mac64.tar.gz')
    elif OS_ID == 'linux':
        if DIST_ID in DEBIAN_DISTS:
            url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/'
                   'sratoolkit.current-ubuntu64.tar.gz')
        elif DIST_ID in REDHAT_DISTS:
            url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/'
                   'sratoolkit.current-centos_linux64.tar.gz')

    dnld_path = os.path.join(DIR_DEP, 'sra-toolkit.tar.gz')

    try:
        fasterq_dump = 'fasterq-dump'
        call(fasterq_dump)
    except Exception:
        try:
            dir_bin = os.path.join(DIR_DEP,
                                   sra_toolkit_dir_name(path=DIR_DEP))
            fasterq_dump = os.path.join(dir_bin, 'bin', 'fasterq-dump')
            call(fasterq_dump)
        except Exception:
            print('\tSRA Toolkit was not found on this system, trying to '
                  'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()

            dir_bin = os.path.join(DIR_DEP,
                                   sra_toolkit_dir_name(path=DIR_DEP))
            fasterq_dump = os.path.join(dir_bin, 'bin', 'fasterq-dump')

            if not os.path.exists(fasterq_dump):
                print('\tCould not download SRA Toolkit.')
                sys.exit(1)

    print('\tSRA Toolkit is available:')
    print('\t\t' + fasterq_dump + '\n')

    return fasterq_dump


# BLAST+
def dep_check_blast(): # noqa
    if OS_ID == 'mac':
        url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/'
               'ncbi-blast-2.9.0+-x64-macosx.tar.gz')
    elif OS_ID == 'linux':
        url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/'
               'ncbi-blast-2.9.0+-x64-linux.tar.gz')

    dnld_path = os.path.join(DIR_DEP, 'ncbi-blast.tar.gz')

    try:
        makeblastdb = 'makeblastdb'
        blastn = 'blastn'
        tblastn = 'tblastn'
        call(makeblastdb)
    except Exception:
        try:
            dir_bin = os.path.join(DIR_DEP, blast_dir_name(path=DIR_DEP))
            makeblastdb = os.path.join(dir_bin, 'bin', 'makeblastdb')
            blastn = os.path.join(dir_bin, 'bin', 'blastn')
            tblastn = os.path.join(dir_bin, 'bin', 'tblastn')
            call(makeblastdb)
        except Exception:
            print('\tBLAST+ was not found on this system, trying to download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()

            dir_bin = os.path.join(DIR_DEP, blast_dir_name(path=DIR_DEP))
            makeblastdb = os.path.join(dir_bin, 'bin', 'makeblastdb')
            blastn = os.path.join(dir_bin, 'bin', 'blastn')
            tblastn = os.path.join(dir_bin, 'bin', 'tblastn')

            if not os.path.exists(makeblastdb) or \
               not os.path.exists(blastn) or \
               not os.path.exists(tblastn):
                print('\tCould not download BLAST+.')
                sys.exit(1)

    print('\tBLAST+ is available:')
    print('\t\t' + makeblastdb)
    print('\t\t' + blastn)
    print('\t\t' + tblastn + '\n')

    return (makeblastdb, blastn, tblastn)


# VSEARCH
def dep_check_vsearch(): # noqa
    url = 'https://github.com/torognes/vsearch/archive/master.tar.gz'
    dnld_path = os.path.join(DIR_DEP, 'vsearch.tar.gz')
    dir_bin = os.path.join(DIR_DEP, 'vsearch-master')

    try:
        vsearch = 'vsearch'
        call(vsearch)
    except Exception:
        try:
            vsearch = os.path.join(dir_bin, 'bin', 'vsearch')
            call(vsearch)
        except Exception:
            print('\tVsearch was not found on this system, trying to '
                  'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()
            try:
                print('\tCompiling Vsearch.')
                call('./autogen.sh', cwd=dir_bin)
                call('./configure', cwd=dir_bin)
                call('make', cwd=dir_bin)
                vsearch = os.path.join(dir_bin, 'bin', 'vsearch')
                call(vsearch)
            except Exception:
                print('\t\tSomething went wrong while trying to compile '
                      'Vsearch.')
                print('\t\tTry downloading and installing it manually from:')
                print('\t\thttps://github.com/torognes/vsearch')
                sys.exit(1)

    print('\tVsearch is available:')
    print('\t\t' + vsearch + '\n')

    return vsearch


# SPAdes
def dep_check_spades(): # noqa
    if OS_ID == 'mac':
        url = ('http://cab.spbu.ru/files/release3.13.1/'
               'SPAdes-3.13.1-Darwin.tar.gz')
    elif OS_ID == 'linux':
        url = ('http://cab.spbu.ru/files/release3.13.1/'
               'SPAdes-3.13.1-Linux.tar.gz')

    dnld_path = os.path.join(DIR_DEP, 'SPAdes.tar.gz')

    try:
        spades = 'spades.py'
        call(spades)
    except Exception:
        try:
            dir_bin = os.path.join(DIR_DEP, spades_dir_name(path=DIR_DEP))
            spades = os.path.join(dir_bin, 'bin', 'spades.py')
            call(spades)
        except Exception:
            print('\tSPAdes was not found on this system, trying to download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()
            try:
                dir_bin = os.path.join(DIR_DEP, spades_dir_name(path=DIR_DEP))
                spades = os.path.join(dir_bin, 'bin', 'spades.py')
                call(spades)
            except Exception:
                print('\tCould not download SPAdes.')
                sys.exit(1)

    print('\tSPAdes is available:')
    print('\t\t' + spades + '\n')

    return spades
