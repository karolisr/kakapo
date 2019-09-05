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
import stat
import tarfile
import zipfile
import re

from kakapo.config import DIR_DEP, DIR_CFG, OS_ID, DIST_ID
from kakapo.helpers import list_of_dirs
from kakapo.helpers import download_file
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


def bowtie2_dir_name(path): # noqa
    ld = list_of_dirs(path=path)
    pattern = 'bowtie2'
    dl = [d for d in ld if pattern in d]
    if len(dl) > 0:
        return dl[0]
    else:
        return ''


def kraken2_dir_name(path): # noqa
    ld = list_of_dirs(path=path)
    pattern = 'kraken2'
    dl = [d for d in ld if pattern in d]
    if len(dl) > 0:
        return dl[0]
    else:
        return ''


# Seqtk
def dep_check_seqtk(logger=print): # noqa
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
            logger('Seqtk was not found on this system, trying to download.')
            download_file(url, dnld_path)
            zip_ref = zipfile.ZipFile(dnld_path, 'r')
            zip_ref.extractall(DIR_DEP)
            zip_ref.close()
            try:
                logger('Compiling Seqtk.')
                call('make', cwd=dir_bin)
                seqtk = os.path.join(dir_bin, 'seqtk')
                call(seqtk)
            except Exception:
                logger('Something went wrong while trying to compile Seqtk.')
                logger('Try downloading and installing it manually from: '
                       'https://github.com/lh3/seqtk')
                sys.exit(1)

    logger('Seqtk is available: ' + seqtk)

    return seqtk


def get_version_seqtk(seqtk):  # noqa
    _, err = call(seqtk)
    v = re.findall(r'Version\:\s([\d\w\.\-]*)', err.decode(),
                   flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


def _write_trimmomatic_adapters_file(logger=print):

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
        logger('Writing Trimmomatic adapter files: ' + path_adapters)
        with open(path_adapters, mode='w') as f:
            f.write(adapters)

    return path_adapters

# Trimmomatic
def dep_check_trimmomatic(logger=print): # noqa
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
        logger('Could not download Trimmomatic.')
        sys.exit(1)

    logger('Trimmomatic is available: ' + trimmomatic)
    path_adapters = _write_trimmomatic_adapters_file()

    return trimmomatic, path_adapters


def get_version_trimmomatic(trimmomatic):  # noqa
    cmd = ['java', '-jar', trimmomatic, '-version']
    out, _ = call(cmd)
    v = out.strip().decode()
    return v


# SRA Toolkit
def dep_check_sra_toolkit(logger=print): # noqa
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
            logger('SRA Toolkit was not found on this system, trying to '
                   'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()

            dir_bin = os.path.join(DIR_DEP,
                                   sra_toolkit_dir_name(path=DIR_DEP))
            fasterq_dump = os.path.join(dir_bin, 'bin', 'fasterq-dump')

            if not os.path.exists(fasterq_dump):
                logger('Could not download SRA Toolkit.')
                sys.exit(1)

    logger('SRA Toolkit is available: ' + fasterq_dump)

    return fasterq_dump


def get_version_fasterq_dump(fasterq_dump):  # noqa
    out, _ = call([fasterq_dump, '--version'])
    v = re.findall(r'\:\s([\d\.]*)', out.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


# BLAST+
def dep_check_blast(logger=print): # noqa
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
            logger('BLAST+ was not found on this system, trying to download.')
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
                logger('Could not download BLAST+.')
                sys.exit(1)

    logger('makeblastdb is available: ' + makeblastdb)
    logger('blastn is available: ' + blastn)
    logger('tblastn is available: ' + tblastn)

    return makeblastdb, blastn, tblastn


def get_version_blast(any_blast_bin):  # noqa
    out, _ = call([any_blast_bin, '-version'])
    v = re.findall(r'\sblast\s([\d\.]*)', out.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


# VSEARCH
def dep_check_vsearch(logger=print): # noqa
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
            logger('Vsearch was not found on this system, trying to download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()
            try:
                logger('Compiling Vsearch.')
                call('./autogen.sh', cwd=dir_bin)
                call('./configure', cwd=dir_bin)
                call('make', cwd=dir_bin)
                vsearch = os.path.join(dir_bin, 'bin', 'vsearch')
                call(vsearch)
            except Exception:
                logger('Something went wrong while trying to compile Vsearch.')
                logger('Try downloading and installing it manually from: '
                       'https://github.com/torognes/vsearch')
                sys.exit(1)

    logger('Vsearch is available: ' + vsearch)

    return vsearch


def get_version_vsearch(vsearch):  # noqa
    _, err = call([vsearch, '-version'])
    v = re.findall(r'^vsearch\sv([\d\.]*)', err.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


# SPAdes
def dep_check_spades(logger=print): # noqa
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
            logger('SPAdes was not found on this system, trying to download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()
            try:
                dir_bin = os.path.join(DIR_DEP, spades_dir_name(path=DIR_DEP))
                spades = os.path.join(dir_bin, 'bin', 'spades.py')
                call(spades)
            except Exception:
                logger('Could not download SPAdes.')
                sys.exit(1)

    logger('SPAdes is available: ' + spades)

    return spades


def get_version_spades(spades):  # noqa
    out, _ = call([spades, '--version'])
    v = re.findall(r'^SPAdes\sv([\d\.]*)', out.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


# Bowtie 2
def dep_check_bowtie2(logger=print): # noqa
    if OS_ID == 'mac':
        url = ('https://sourceforge.net/projects/bowtie-bio/files/bowtie2/'
               '2.3.5.1/bowtie2-2.3.5.1-macos-x86_64.zip/download')
    elif OS_ID == 'linux':
        url = ('https://sourceforge.net/projects/bowtie-bio/files/bowtie2/'
               '2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip/download')

    dnld_path = os.path.join(DIR_DEP, 'bowtie2.zip')

    try:
        bowtie2 = 'bowtie2'
        bowtie2_build = 'bowtie2-build'
        call(bowtie2)
        call(bowtie2_build)
    except Exception:
        try:
            dir_bin = os.path.join(DIR_DEP, bowtie2_dir_name(path=DIR_DEP))
            bowtie2 = os.path.join(dir_bin, 'bowtie2')
            bowtie2_build = os.path.join(dir_bin, 'bowtie2-build')
            call(bowtie2)
            call(bowtie2_build)
        except Exception:
            logger('Bowtie 2 was not found on this system, trying to '
                   'download.')
            download_file(url, dnld_path)
            zip_ref = zipfile.ZipFile(dnld_path, 'r')
            zip_ref.extractall(DIR_DEP)
            zip_ref.close()

            dir_bin = os.path.join(DIR_DEP, bowtie2_dir_name(path=DIR_DEP))
            bowtie2 = os.path.join(dir_bin, 'bowtie2')
            bowtie2_build = os.path.join(dir_bin, 'bowtie2-build')

            bowtie2_execs = ('',
                             '-align-l',
                             '-align-l-debug',
                             '-align-s',
                             '-align-s-debug',
                             '-build',
                             '-build-l',
                             '-build-l-debug',
                             '-build-s',
                             '-build-s-debug',
                             '-inspect',
                             '-inspect-l',
                             '-inspect-l-debug',
                             '-inspect-s',
                             '-inspect-s-debug')

            for bt2exe in bowtie2_execs:
                os.chmod(bowtie2 + bt2exe, stat.S_IRWXU | stat.S_IRGRP |
                         stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)

            if not os.path.exists(bowtie2):
                logger('Could not download Bowtie 2.')
                sys.exit(1)

    logger('Bowtie 2 is available: ' + bowtie2)

    return bowtie2, bowtie2_build


def get_version_bowtie2(bowtie2):  # noqa
    out, _ = call([bowtie2, '--version'])
    v = re.findall(r'^.*?version\s([\d\.]*)', out.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


# Kraken 2
def dep_check_kraken2(logger=print): # noqa
    url = ('https://github.com/DerrickWood/kraken2/archive/v2.0.8-beta.tar.gz')

    dnld_path = os.path.join(DIR_DEP, 'kraken2.tar.gz')

    try:
        kraken2 = 'kraken2'
        kraken2_build = 'kraken2-build'
        call(kraken2)
        call(kraken2_build)
    except Exception:
        try:
            dir_bin = os.path.join(DIR_DEP, kraken2_dir_name(path=DIR_DEP))
            kraken2 = os.path.join(dir_bin, 'kraken2')
            kraken2_build = os.path.join(dir_bin, 'kraken2-build')
            call(kraken2)
            call(kraken2_build)
        except Exception:
            logger('Kraken 2 was not found on this system, trying to '
                   'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()

            dir_bin = os.path.join(DIR_DEP, kraken2_dir_name(path=DIR_DEP))
            kraken2 = os.path.join(dir_bin, 'kraken2')
            kraken2_build = os.path.join(dir_bin, 'kraken2-build')

            try:
                logger('Compiling Kraken 2')
                call(['./install_kraken2.sh', '.'], cwd=dir_bin)
            except Exception:
                logger('Something went wrong while trying to compile'
                       'Kraken 2.')
                logger('Try downloading and installing it manually from: '
                       'https://github.com/DerrickWood/kraken2')

            if not os.path.exists(kraken2):
                logger('Could not download Kraken 2.')
                sys.exit(1)

    logger('Kraken 2 is available: ' + kraken2)

    return kraken2, kraken2_build


def get_version_kraken2(kraken2):  # noqa
    out, _ = call([kraken2, '--version'])
    v = re.findall(r'^.*?version\s([\d\.\-A-Za-z]*)', out.decode(),
                   flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v
