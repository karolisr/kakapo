# -*- coding: utf-8 -*-

"""
Dependencies.
"""

import os
import re
import stat
import sys
import tarfile
import zipfile

from os.path import join as opj
from shutil import move

from kakapo.config import DIR_DEP, DIR_CFG, OS_ID, DIST_ID
from kakapo.helpers import download_file
from kakapo.helpers import list_of_dirs
from kakapo.helpers import replace_line_in_file
from kakapo.os_diffs import DEBIAN_DISTS, REDHAT_DISTS
from kakapo.shell import call

VER_UNK = '[version could not be determined]'


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


def rcorrector_dir_name(path): # noqa
    ld = list_of_dirs(path=path)
    pattern = 'Rcorrector'
    dl = [d for d in ld if pattern in d]
    if len(dl) > 0:
        return dl[0]
    else:
        return ''


def get_version_seqtk(seqtk):  # noqa
    _, err = call(seqtk)
    if type(err) not in (bytes, str) or err == b'':
        return VER_UNK
    v = re.findall(r'Version\:\s([\d\w\.\-]*)', err.decode(),
                   flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


def get_version_trimmomatic(trimmomatic):  # noqa
    cmd = ['java', '-jar', trimmomatic, '-version']
    out, _ = call(cmd)
    if type(out) not in (bytes, str) or out == b'':
        return VER_UNK
    v = out.strip().decode()
    return v


def get_version_fasterq_dump(fasterq_dump):  # noqa
    out, _ = call([fasterq_dump, '--version'])
    if type(out) not in (bytes, str) or out == b'':
        return VER_UNK
    v = re.findall(r'[\:n]\s([\d\.]*)', out.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


def get_version_blast(any_blast_bin):  # noqa
    out, _ = call([any_blast_bin, '-version'])
    if type(out) not in (bytes, str) or out == b'':
        return VER_UNK
    v = re.findall(r'\sblast\s([\d\.]*)', out.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


def get_version_vsearch(vsearch):  # noqa
    _, err = call([vsearch, '-version'])
    if type(err) not in (bytes, str) or err == b'':
        return VER_UNK
    v = re.findall(r'^vsearch\sv([\d\.]*)', err.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


def get_version_spades(spades):  # noqa
    out, _ = call([spades, '--version'])
    if type(out) not in (bytes, str) or out == b'':
        return VER_UNK
    v = re.findall(r'^SPAdes.*v([\d\.]*)', out.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


def get_version_bowtie2(bowtie2):  # noqa
    out, _ = call([bowtie2, '--version'])
    if type(out) not in (bytes, str) or out == b'':
        return VER_UNK
    v = re.findall(r'^.*?version\s([\d\.]*)', out.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


def get_version_rcorrector(rcorrector):  # noqa
    out, _ = call([rcorrector, '-version'])
    if type(out) not in (bytes, str) or out == b'':
        return VER_UNK
    v = re.findall(r'^Rcorrector\sv([\d\.]*)', out.decode(), flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


def get_version_kraken2(kraken2):  # noqa
    out, _ = call([kraken2, '--version'])
    if type(out) not in (bytes, str) or out == b'':
        return VER_UNK
    v = re.findall(r'^.*?version\s([\d\.\-A-Za-z]*)', out.decode(),
                   flags=re.MULTILINE)
    if len(v) > 0:
        v = v[0]
    return v


# Seqtk
def dep_check_seqtk(force=False, logger=print): # noqa
    url = 'https://github.com/lh3/seqtk/archive/master.zip'
    dnld_path = opj(DIR_DEP, 'seqtk.zip')
    dir_bin = opj(DIR_DEP, 'seqtk-master')

    try:
        if force is True:
            raise
        seqtk = 'seqtk'
        call(seqtk)
    except Exception:
        try:
            seqtk = opj(dir_bin, 'seqtk')
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
                seqtk = opj(dir_bin, 'seqtk')
                call(seqtk)
            except Exception:
                replace_line_in_file(opj(dir_bin, 'Makefile'), 'CC=gcc', 'CC=cc')
                try:
                    call('make', cwd=dir_bin)
                    seqtk = opj(dir_bin, 'seqtk')
                    call(seqtk)
                except Exception:
                    logger('Something went wrong while trying to compile Seqtk.')
                    logger('Try downloading and installing it manually from: '
                           'https://github.com/lh3/seqtk')
                    return None

    v = get_version_seqtk(seqtk)
    logger('Seqtk is available: ' + v + ' ' + seqtk)

    return seqtk


def _write_trimmomatic_adapters_file(logger=print): # noqa

    path_adapters = opj(DIR_DEP, 'trimmomatic_adapters.fasta')

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
    dnld_path = opj(DIR_DEP, 'Trimmomatic-0.39.zip')
    dir_bin = opj(DIR_DEP, 'Trimmomatic-0.39')
    trimmomatic = opj(dir_bin, 'trimmomatic-0.39.jar')

    if not os.path.exists(trimmomatic):
        download_file(url, dnld_path)
        zip_ref = zipfile.ZipFile(dnld_path, 'r')
        zip_ref.extractall(DIR_DEP)
        zip_ref.close()

    if not os.path.exists(trimmomatic):
        logger('Could not download Trimmomatic.')
        return None, None

    v = get_version_trimmomatic(trimmomatic)
    logger('Trimmomatic is available: ' + v + ' ' + trimmomatic)

    path_adapters = _write_trimmomatic_adapters_file(logger)

    return trimmomatic, path_adapters


# SRA Toolkit
def dep_check_sra_toolkit(force=False, logger=print): # noqa
    if OS_ID == 'mac':
        url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/'
               'sratoolkit.2.9.6-mac64.tar.gz')
    elif OS_ID == 'linux':
        if DIST_ID in DEBIAN_DISTS:
            url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.2/'
                   'sratoolkit.2.10.2-ubuntu64.tar.gz')
        elif DIST_ID in REDHAT_DISTS:
            url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.2/'
                   'sratoolkit.2.10.2-centos_linux64.tar.gz')

    dnld_path = opj(DIR_DEP, 'sra-toolkit.tar.gz')

    try:
        if force is True:
            raise
        fasterq_dump = 'fasterq-dump'
        call(fasterq_dump)
    except Exception:
        try:
            dir_bin = opj(DIR_DEP, sra_toolkit_dir_name(path=DIR_DEP))
            fasterq_dump = opj(dir_bin, 'bin', 'fasterq-dump')
            call(fasterq_dump)
        except Exception:
            logger('SRA Toolkit was not found on this system, trying to '
                   'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()

            dir_bin = opj(DIR_DEP, sra_toolkit_dir_name(path=DIR_DEP))
            fasterq_dump = opj(dir_bin, 'bin', 'fasterq-dump')

            if not os.path.exists(fasterq_dump):
                logger('Could not download SRA Toolkit.')
                return None

    v = get_version_fasterq_dump(fasterq_dump)
    logger('fasterq-dump is available: ' + v + ' ' + fasterq_dump)

    return fasterq_dump


# BLAST+
def dep_check_blast(force=False, logger=print): # noqa
    if OS_ID == 'mac':
        url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/'
               'ncbi-blast-2.10.0+-x64-macosx.tar.gz')
    elif OS_ID == 'linux':
        if DIST_ID in DEBIAN_DISTS:
            url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/'
                   'ncbi-blast-2.10.0+-x64-linux.tar.gz')
        elif DIST_ID in REDHAT_DISTS:
            url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/'
                   'ncbi-blast-2.10.0+-x64-linux.tar.gz')

    dnld_path = opj(DIR_DEP, 'ncbi-blast.tar.gz')

    try:
        if force is True:
            raise
        makeblastdb = 'makeblastdb'
        blastn = 'blastn'
        tblastn = 'tblastn'
        call(makeblastdb)
    except Exception:
        try:
            dir_bin = opj(DIR_DEP, blast_dir_name(path=DIR_DEP))
            makeblastdb = opj(dir_bin, 'bin', 'makeblastdb')
            blastn = opj(dir_bin, 'bin', 'blastn')
            tblastn = opj(dir_bin, 'bin', 'tblastn')
            call(makeblastdb)
        except Exception:
            logger('BLAST+ was not found on this system, trying to download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()

            dir_bin = opj(DIR_DEP, blast_dir_name(path=DIR_DEP))
            makeblastdb = opj(dir_bin, 'bin', 'makeblastdb')
            blastn = opj(dir_bin, 'bin', 'blastn')
            tblastn = opj(dir_bin, 'bin', 'tblastn')

            if not os.path.exists(makeblastdb) or \
               not os.path.exists(blastn) or \
               not os.path.exists(tblastn):
                logger('Could not download BLAST+.')
                return None, None, None

    v = get_version_blast(makeblastdb)
    logger('makeblastdb is available: ' + v + ' ' + makeblastdb)
    v = get_version_blast(blastn)
    logger('blastn is available: ' + v + ' ' + blastn)
    v = get_version_blast(tblastn)
    logger('tblastn is available: ' + v + ' ' + tblastn)

    return makeblastdb, blastn, tblastn


# VSEARCH
def dep_check_vsearch(force=False, logger=print): # noqa
    url = 'https://github.com/torognes/vsearch/archive/master.tar.gz'
    dnld_path = opj(DIR_DEP, 'vsearch.tar.gz')
    dir_bin = opj(DIR_DEP, 'vsearch-master')

    try:
        if force is True:
            raise
        vsearch = 'vsearch'
        call(vsearch)
    except Exception:
        try:
            vsearch = opj(dir_bin, 'bin', 'vsearch')
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
                vsearch = opj(dir_bin, 'bin', 'vsearch')
                call(vsearch)
            except Exception:
                logger('Something went wrong while trying to compile Vsearch.')
                logger('Try downloading and installing it manually from: '
                       'https://github.com/torognes/vsearch')
                return None

    v = get_version_vsearch(vsearch)
    logger('Vsearch is available: ' + v + ' ' + vsearch)

    return vsearch


# SPAdes
def dep_check_spades(force=False, logger=print): # noqa
    if OS_ID == 'mac':
        url = ('http://cab.spbu.ru/files/release3.14.0/'
               'SPAdes-3.14.0-Darwin.tar.gz')
    elif OS_ID == 'linux':
        url = ('http://cab.spbu.ru/files/release3.14.0/'
               'SPAdes-3.14.0-Linux.tar.gz')

    dnld_path = opj(DIR_DEP, 'SPAdes.tar.gz')

    try:
        if force is True:
            raise
        spades = 'spades.py'
        call(spades)
    except Exception:
        try:
            dir_bin = opj(DIR_DEP, spades_dir_name(path=DIR_DEP))
            spades = opj(dir_bin, 'bin', 'spades.py')
            call(spades)
        except Exception:
            logger('SPAdes was not found on this system, trying to download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()
            try:
                dir_bin = opj(DIR_DEP, spades_dir_name(path=DIR_DEP))
                spades = opj(dir_bin, 'bin', 'spades.py')
                call(spades)
            except Exception:
                logger('Could not download SPAdes.')
                return None

    v = get_version_spades(spades)
    logger('SPAdes is available: ' + v + ' ' + spades)

    return spades


# Bowtie 2
def dep_check_bowtie2(force=False, logger=print): # noqa
    if OS_ID == 'mac':
        url = ('https://sourceforge.net/projects/bowtie-bio/files/bowtie2/'
               '2.3.5.1/bowtie2-2.3.5.1-macos-x86_64.zip/download')
    elif OS_ID == 'linux':
        url = ('https://sourceforge.net/projects/bowtie-bio/files/bowtie2/'
               '2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip/download')

    dnld_path = opj(DIR_DEP, 'bowtie2.zip')

    try:
        if force is True:
            raise
        bowtie2 = 'bowtie2'
        bowtie2_build = 'bowtie2-build'
        call(bowtie2)
        call(bowtie2_build)
    except Exception:
        try:
            dir_bin = opj(DIR_DEP, bowtie2_dir_name(path=DIR_DEP))
            bowtie2 = opj(dir_bin, 'bowtie2')
            bowtie2_build = opj(dir_bin, 'bowtie2-build')
            call(bowtie2)
            call(bowtie2_build)
        except Exception:
            logger('Bowtie 2 was not found on this system, trying to '
                   'download.')
            download_file(url, dnld_path)
            zip_ref = zipfile.ZipFile(dnld_path, 'r')
            zip_ref.extractall(DIR_DEP)
            zip_ref.close()

            dir_bin = opj(DIR_DEP, bowtie2_dir_name(path=DIR_DEP))
            bowtie2 = opj(dir_bin, 'bowtie2')
            bowtie2_build = opj(dir_bin, 'bowtie2-build')

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
                return None, None

    v = get_version_bowtie2(bowtie2)
    logger('bowtie2 is available: ' + v + ' ' + bowtie2)
    v = get_version_bowtie2(bowtie2_build)
    logger('bowtie2-build is available: ' + v + ' ' + bowtie2_build)

    return bowtie2, bowtie2_build


# Rcorrector
def dep_check_rcorrector(force=False, logger=print): # noqa
    url = 'https://github.com/mourisl/Rcorrector/archive/master.tar.gz'
    dnld_path = opj(DIR_DEP, 'rcorrector.tar.gz')

    try:
        if force is True:
            raise
        rcorrector = 'run_rcorrector.pl'
        call(rcorrector)
    except Exception:
        try:
            dir_bin = opj(DIR_DEP, rcorrector_dir_name(path=DIR_DEP))
            rcorrector = opj(dir_bin, 'run_rcorrector.pl')
            call(rcorrector)
        except Exception:
            logger('Rcorrector was not found on this system, trying to '
                   'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()
            dir_bin = opj(DIR_DEP, rcorrector_dir_name(path=DIR_DEP))
            try:
                logger('Compiling Rcorrector.')
                call('make', cwd=dir_bin)
                rcorrector = opj(dir_bin, 'run_rcorrector.pl')
                os.chmod(rcorrector, stat.S_IRWXU | stat.S_IRGRP |
                         stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
                call(rcorrector)
            except Exception:
                logger('Something went wrong while trying to compile '
                       'Rcorrector.')
                logger('Try downloading and installing it manually from: '
                       'https://github.com/mourisl/Rcorrector')
                return None

    v = get_version_rcorrector(rcorrector)
    logger('Rcorrector is available: ' + v + ' ' + rcorrector)

    return rcorrector


# Kraken2
def dep_check_kraken2(force=False, logger=print): # noqa
    url = ('https://github.com/karolisr/kraken2/archive/master.tar.gz')

    dnld_path = opj(DIR_DEP, 'kraken2.tar.gz')

    try:
        if force is True:
            raise
        kraken2 = 'kraken2'
        kraken2_build = 'kraken2-build'
        call(kraken2)
        call(kraken2_build)
    except Exception:
        try:
            dir_bin = opj(DIR_DEP, kraken2_dir_name(path=DIR_DEP))
            kraken2 = opj(dir_bin, 'bin', 'kraken2')
            kraken2_build = opj(dir_bin, 'bin', 'kraken2-build')
            call(kraken2)
            call(kraken2_build)
        except Exception:
            logger('Kraken2 was not found on this system, trying to '
                   'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(DIR_DEP)
            tar_ref.close()

            dir_bin = opj(DIR_DEP, kraken2_dir_name(path=DIR_DEP))
            kraken2 = opj(dir_bin, 'bin', 'kraken2')
            kraken2_build = opj(dir_bin, 'bin', 'kraken2-build')

            try:
                logger('Compiling Kraken2')
                call(['./install_kraken2.sh', 'bin'], cwd=dir_bin)
            except Exception:
                pass

            if not os.path.exists(kraken2):
                logger('Something went wrong while trying to compile '
                       'Kraken2.')
                logger('Try downloading and installing it manually from: '
                       'https://github.com/DerrickWood/kraken2')
                return None, None

    v = get_version_kraken2(kraken2)
    logger('kraken2 is available: ' + v + ' ' + kraken2)
    v = get_version_kraken2(kraken2_build)
    logger('kraken2-build is available: ' + v + ' ' + kraken2_build)

    return kraken2, kraken2_build


def download_kraken2_dbs(dbs_path, dnld):  # noqa

    base_kraken2_url = 'ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/'

    if dnld is True:

        # --------------------------------------------------------------------

        base = '16S_Silva_20190418'
        url = base_kraken2_url + base + '.tgz'
        tgz = opj(dbs_path, base + '.tgz')
        p_orig = opj(dbs_path, base)
        db_name = '16S_Silva'
        p_new = opj(dbs_path, db_name)

        if not os.path.exists(p_new):
            download_file(url=url, local_path=tgz, protocol='ftp')
            tar_ref = tarfile.open(tgz, 'r:gz')
            tar_ref.extractall(dbs_path)
            tar_ref.close()
            os.remove(tgz)
            move(p_orig, p_new)

        # --------------------------------------------------------------------

        base = 'minikraken2_v2_8GB_201904_UPDATE'
        url = base_kraken2_url + base + '.tgz'
        tgz = opj(dbs_path, base + '.tgz')
        p_orig = opj(dbs_path, base)
        db_name = 'minikraken2_v2'
        p_new = opj(dbs_path, db_name)

        if not os.path.exists(p_new):
            download_file(url=url, local_path=tgz, protocol='ftp')
            tar_ref = tarfile.open(tgz, 'r:gz')
            tar_ref.extractall(dbs_path)
            tar_ref.close()
            os.remove(tgz)
            move(p_orig, p_new)

        # --------------------------------------------------------------------

        base_dropbox_url = 'https://www.dropbox.com/s/'

        base = 'mitochondrion_and_plastid'
        garb = 'vkbp7iys6s76tvf/'
        url = base_dropbox_url + garb + base + '.tar.gz?dl=1'
        tgz = opj(dbs_path, base + '.tar.gz')
        p = opj(dbs_path, base)

        if not os.path.exists(p):
            download_file(url=url, local_path=tgz, protocol='http')
            tar_ref = tarfile.open(tgz, 'r:gz')
            tar_ref.extractall(dbs_path)
            tar_ref.close()
            os.remove(tgz)

        # --------------------------------------------------------------------

        base_dropbox_url = 'https://www.dropbox.com/s/'

        base = 'mitochondrion'
        garb = '6liwneb26uvjuec/'
        url = base_dropbox_url + garb + base + '.tar.gz?dl=1'
        tgz = opj(dbs_path, base + '.tar.gz')
        p = opj(dbs_path, base)

        if not os.path.exists(p):
            download_file(url=url, local_path=tgz, protocol='http')
            tar_ref = tarfile.open(tgz, 'r:gz')
            tar_ref.extractall(dbs_path)
            tar_ref.close()
            os.remove(tgz)

        # --------------------------------------------------------------------

        base_dropbox_url = 'https://www.dropbox.com/s/'

        base = 'plastid'
        garb = 's9vdg4mxrfy1szn/'
        url = base_dropbox_url + garb + base + '.tar.gz?dl=1'
        tgz = opj(dbs_path, base + '.tar.gz')
        p = opj(dbs_path, base)

        if not os.path.exists(p):
            download_file(url=url, local_path=tgz, protocol='http')
            tar_ref = tarfile.open(tgz, 'r:gz')
            tar_ref.extractall(dbs_path)
            tar_ref.close()
            os.remove(tgz)

        # --------------------------------------------------------------------

        base_dropbox_url = 'https://www.dropbox.com/s/'

        base = 'viral'
        garb = '7xz31c7vw088n27/'
        url = base_dropbox_url + garb + base + '.tar.gz?dl=1'
        tgz = opj(dbs_path, base + '.tar.gz')
        p = opj(dbs_path, base)

        if not os.path.exists(p):
            download_file(url=url, local_path=tgz, protocol='http')
            tar_ref = tarfile.open(tgz, 'r:gz')
            tar_ref.extractall(dbs_path)
            tar_ref.close()
            os.remove(tgz)

        # --------------------------------------------------------------------

    kraken2_dbs = {k: opj(dbs_path, k) for k in list_of_dirs(dbs_path)}
    return kraken2_dbs
