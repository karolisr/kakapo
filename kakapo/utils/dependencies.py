"""Install Kakapo dependencies."""

import stat
import tarfile
import zipfile

from os import chmod
from os import linesep as lns
from os import remove
from os.path import basename
from os.path import dirname
from os.path import exists as ope
from os.path import join as opj
from shutil import move
from shutil import rmtree
from tempfile import NamedTemporaryFile

from kakapo.utils.c.kakapolib import dep_check_kakapolib as dckkpl
from kakapo.utils.homebrew import brew_get
from kakapo.utils.logging import Log
from kakapo.utils.misc import download_file
from kakapo.utils.misc import list_of_dirs_at_path
from kakapo.utils.misc import replace_line_in_file
from kakapo.utils.subp import run
from kakapo.utils.subp import run_then_grep
from kakapo.utils.subp import which


PY3 = which('python3')


def dep_check_kakapolib(force):
    return dckkpl(force)


def dnld_kraken2_dbs(dbs_path):
    Log.inf('Checking for available Kraken2 databases.')
    kraken2_dbs = download_kraken2_dbs(dbs_path)
    for db in sorted(kraken2_dbs.keys()):
        Log.msg('Found Kraken2 database:', db)
    return kraken2_dbs


def get_dep_version(cmd, regexp, do_not_raise=True):
    v = run_then_grep(cmd=cmd, regexp=regexp, do_not_raise=do_not_raise)
    if len(v) > 0:
        return v[0]
    else:
        return '?'


def get_dep_dir(path, pattern):
    ld, e = list_of_dirs_at_path(path)
    if ld is not None:
        dd = [d for d in ld if pattern in d]
        if len(dd) > 0:
            return dd[0]
        else:
            return ''


# Seqtk
def dep_check_seqtk(dir_dep, force):
    url = 'https://github.com/lh3/seqtk/archive/master.zip'
    dnld_path = opj(dir_dep, 'seqtk.zip')
    dir_bin = opj(dir_dep, 'seqtk-master')

    fp = NamedTemporaryFile()
    fp.write(str.encode('>seq' + lns + 'ATGC'))
    fp.seek(0)
    cmd = ['', 'seq', '-r', fp.name]

    try:
        if force is True:
            raise
        seqtk = which('seqtk')
        cmd[0] = seqtk
        run(cmd, do_not_raise=True)
    except Exception:
        try:
            seqtk = opj(dir_bin, 'seqtk')
            cmd[0] = seqtk
            run(cmd, do_not_raise=True)
        except Exception:
            Log.wrn('Seqtk was not found on this system, trying to download.')
            download_file(url, dnld_path)
            zip_ref = zipfile.ZipFile(dnld_path, 'r')
            zip_ref.extractall(dir_dep)
            zip_ref.close()
            try:
                Log.wrn('Compiling Seqtk.')
                run('make', cwd=dir_bin)
                run(cmd, do_not_raise=True)
            except Exception:
                replace_line_in_file(opj(dir_bin, 'Makefile'),
                                     'CC=gcc', 'CC=cc')
                try:
                    run('make', cwd=dir_bin)
                    run(cmd, do_not_raise=True)
                except Exception:
                    Log.err(
                        'Something went wrong while trying to compile Seqtk.')
                    Log.msg('Try downloading and installing it manually from: '
                            'https://github.com/lh3/seqtk')
                    fp.close()
                    return None

    fp.close()

    v = get_dep_version([seqtk], r'Version\:\s([\d\w\.\-]*)')
    Log.msg('Seqtk is available:', v + ' ' + seqtk)

    return seqtk


def _write_trimmomatic_adapters_file(dir_dep):
    path_adapters = opj(dir_dep, 'trimmomatic_adapters.fasta')

    adapters = ('>TruSeq2_SE\n'
                'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG\n'
                '>TruSeq2_PE_f\n'
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n'
                '>TruSeq2_PE_r\n'
                'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG\n'
                '>TruSeq3_IndexedAdapter\n'
                'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n'
                '>TruSeq3_UniversalAdapter\n'
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n'
                '>PrefixPE/1\n'
                'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n'
                '>PrefixPE/2\n'
                'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT\n'
                '>PCR_Primer1\n'
                'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n'
                '>PCR_Primer1_rc\n'
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\n'
                '>PCR_Primer2\n'
                'CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT\n'
                '>PCR_Primer2_rc\n'
                'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG\n'
                '>FlowCell1\n'
                'TTTTTTTTTTAATGATACGGCGACCACCGAGATCTACAC\n'
                '>FlowCell2\n'
                'TTTTTTTTTTCAAGCAGAAGACGGCATACGA\n'
                '>PrefixPE/1\n'
                'TACACTCTTTCCCTACACGACGCTCTTCCGATCT\n'
                '>PrefixPE/2\n'
                'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n'
                '>PE1\n'
                'TACACTCTTTCCCTACACGACGCTCTTCCGATCT\n'
                '>PE1_rc\n'
                'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA\n'
                '>PE2\n'
                'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n'
                '>PE2_rc\n'
                'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n')

    if not ope(path_adapters):
        Log.msg('Writing Trimmomatic adapter files: ' + path_adapters)
        with open(path_adapters, mode='w') as f:
            f.write(adapters)

    return path_adapters


# Trimmomatic
def dep_check_trimmomatic(dir_dep):
    url = ('http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/'
           'Trimmomatic-0.39.zip')
    dnld_path = opj(dir_dep, 'Trimmomatic-0.39.zip')
    dir_bin = opj(dir_dep, 'Trimmomatic-0.39')
    trimmomatic = opj(dir_bin, 'trimmomatic-0.39.jar')

    if not ope(trimmomatic):
        download_file(url, dnld_path)
        zip_ref = zipfile.ZipFile(dnld_path, 'r')
        zip_ref.extractall(dir_dep)
        zip_ref.close()

    if not ope(trimmomatic):
        Log.err('Could not download Trimmomatic.')
        return None, None

    v = get_dep_version(['java', '-jar', trimmomatic, '-version'],
                        r'\d+\.\d+')
    Log.msg('Trimmomatic is available:', v + ' ' + trimmomatic)

    path_adapters = _write_trimmomatic_adapters_file(dir_dep)

    return trimmomatic, path_adapters


# SRA Toolkit
def _ensure_vdb_cfg(dir_bin):
    """
    Ensure that the required configuration files are created without user
    interaction 'vdb-config --interactive'.

    Solves this problem:
        This sra toolkit installation has not been configured.
        Before continuing, please run: vdb-config --interactive
        For more information, see https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud/

    """
    vdb_config = opj(dir_bin, 'bin', 'vdb-config')
    run([vdb_config, '--interactive'], in_txt='x', do_not_raise=True)


def dep_check_sra_toolkit(dir_dep, os_id, dist_id, debian_dists, redhat_dists,
                          force):
    if os_id == 'mac':
        url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/'
               'sratoolkit.2.11.3-mac64.tar.gz')
    elif os_id == 'linux':
        if dist_id in debian_dists:
            url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/'
                   'sratoolkit.2.11.3-ubuntu64.tar.gz')
        elif dist_id in redhat_dists:
            url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/'
                   'sratoolkit.2.11.3-centos_linux64.tar.gz')
        else:
            url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/'
                   'sratoolkit.2.11.3-ubuntu64.tar.gz')

    dnld_path = opj(dir_dep, 'sra-toolkit.tar.gz')

    fasterq_dump = None
    try:
        if force is True:
            raise
        fasterq_dump = which('fasterq-dump')
        dir_bin = dirname(fasterq_dump).strip('bin')
        _ensure_vdb_cfg(dir_bin)
        run([fasterq_dump, "-V"])
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'sratoolkit'))
            _ensure_vdb_cfg(dir_bin)
            fasterq_dump = opj(dir_bin, 'bin', 'fasterq-dump')
            run([fasterq_dump, "-V"])
        except Exception:
            Log.wrn('SRA Toolkit was not found on this system, trying to '
                    'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(dir_dep)
            tar_ref.close()

            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'sratoolkit'))
            fasterq_dump = opj(dir_bin, 'bin', 'fasterq-dump')

            _ensure_vdb_cfg(dir_bin)

            if not ope(fasterq_dump):
                Log.err('Could not download SRA Toolkit.')
                return None

    v = get_dep_version([fasterq_dump, '--version'], r':\s([\d\.]*)')
    if v == '?':
        v = get_dep_version([fasterq_dump, '--version'], r'version\s([\d\.]*)')
    Log.msg('fasterq-dump is available:', v + ' ' + fasterq_dump)

    return fasterq_dump


# BLAST+
def dep_check_blast(dir_dep, os_id, dist_id, debian_dists, redhat_dists,
                    force):
    if os_id == 'mac':
        url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/'
               '2.12.0/ncbi-blast-2.10.1+-x64-macosx.tar.gz')
    elif os_id == 'linux':
        if dist_id in debian_dists:
            url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/'
                   '2.12.0/ncbi-blast-2.12.0+-x64-linux.tar.gz')
        elif dist_id in redhat_dists:
            url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/'
                   '2.12.0/ncbi-blast-2.12.0+-x64-linux.tar.gz')
        else:
            url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/'
                   '2.12.0/ncbi-blast-2.12.0+-x64-linux.tar.gz')


    dnld_path = opj(dir_dep, 'ncbi-blast.tar.gz')

    makeblastdb = None
    blastn = None
    tblastn = None

    try:
        if force is True:
            raise
        makeblastdb = which('makeblastdb')
        blastn = which('blastn')
        tblastn = which('tblastn')
        run([makeblastdb, '-help'])
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'ncbi-blast'))
            makeblastdb = opj(dir_bin, 'bin', 'makeblastdb')
            blastn = opj(dir_bin, 'bin', 'blastn')
            tblastn = opj(dir_bin, 'bin', 'tblastn')
            run([makeblastdb, '-help'])
        except Exception:
            Log.wrn('BLAST+ was not found on this system, trying to download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(dir_dep)
            tar_ref.close()

            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'ncbi-blast'))
            makeblastdb = opj(dir_bin, 'bin', 'makeblastdb')
            blastn = opj(dir_bin, 'bin', 'blastn')
            tblastn = opj(dir_bin, 'bin', 'tblastn')

            if not ope(makeblastdb) or \
                    not ope(blastn) or \
                    not ope(tblastn):
                Log.err('Could not download BLAST+.')
                return None, None, None

    regexp = r'\sblast\s([\d\.]*)'
    v = get_dep_version([makeblastdb, '-version'], regexp)
    Log.msg('makeblastdb is available:', v + ' ' + makeblastdb)
    v = get_dep_version([blastn, '-version'], regexp)
    Log.msg('blastn is available:', v + ' ' + blastn)
    v = get_dep_version([tblastn, '-version'], regexp)
    Log.msg('tblastn is available:', v + ' ' + tblastn)

    return makeblastdb, blastn, tblastn


# VSEARCH
def dep_check_vsearch(dir_dep, os_id, dist_id, debian_dists, redhat_dists,
                      force):
    if os_id == 'mac':
        url = ('https://github.com/torognes/vsearch/releases/download/v2.21.1/'
               'vsearch-2.21.1-macos-x86_64.tar.gz')
    elif os_id == 'linux':
        if dist_id in debian_dists:
            url = ('https://github.com/torognes/vsearch/releases/download/'
                   'v2.21.1/vsearch-2.21.1-linux-x86_64-static.tar.gz')
        elif dist_id in redhat_dists:
            url = ('https://github.com/torognes/vsearch/releases/download/'
                   'v2.21.1/vsearch-2.21.1-linux-x86_64-static.tar.gz')
        else:
            url = ('https://github.com/torognes/vsearch/releases/download/'
                   'v2.21.1/vsearch-2.21.1-linux-x86_64-static.tar.gz')

    dnld_path = opj(dir_dep, 'vsearch.tar.gz')

    try:
        if force is True:
            raise
        vsearch = which('vsearch')
        run(vsearch)
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'vsearch'))
            vsearch = opj(dir_bin, 'bin', 'vsearch')
            run(vsearch)
        except Exception:
            Log.wrn(
                'Vsearch was not found on this system, trying to download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(dir_dep)
            tar_ref.close()
            try:
                dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'vsearch'))
                vsearch = opj(dir_bin, 'bin', 'vsearch')
                if not ope(vsearch):
                    Log.err('Could not download Vsearch.')
                    return None
                else:
                    run(vsearch)
            except Exception:
                Log.err('Vsearch was downloaded, but does not execute.')
                Log.msg('Try downloading and installing it manually from: '
                        'https://github.com/torognes/vsearch')
                return None

    v = get_dep_version([vsearch, '-version'], r'vsearch\sv([\d\.]*)')
    Log.msg('Vsearch is available:', v + ' ' + vsearch)

    return vsearch


# SPAdes
def dep_check_spades(dir_dep, os_id, force):
    if os_id == 'mac':
        url = ('https://github.com/ablab/spades/releases/download/v3.15.4/'
               'SPAdes-3.15.4-Darwin.tar.gz')
    elif os_id == 'linux':
        url = ('https://github.com/ablab/spades/releases/download/v3.15.4/'
               'SPAdes-3.15.4-Linux.tar.gz')

    dnld_path = opj(dir_dep, 'SPAdes.tar.gz')

    try:
        if force is True:
            raise
        spades = which('spades.py')
        run([PY3, spades])
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'SPAdes'))
            spades = opj(dir_bin, 'bin', 'spades.py')
            run([PY3, spades])
        except Exception:
            Log.wrn('SPAdes was not found on this system, trying to download.')
            try:
                download_file(url, dnld_path)
                tar_ref = tarfile.open(dnld_path, 'r:gz')
                tar_ref.extractall(dir_dep)
                tar_ref.close()
            except Exception:
                Log.err('Could not download SPAdes.')
                return None
            try:
                dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'SPAdes'))
                spades = opj(dir_bin, 'bin', 'spades.py')
                # replace_line_in_file(spades,
                #                      '#!/usr/bin/env python',
                #                      '#!/usr/bin/env python3')
                if ope(spades):
                    run([PY3, spades])
                else:
                    Log.err('Could not download SPAdes.')
                    return None
            except Exception:
                Log.err('SPAdes was downloaded, but does not execute.')
                return None

    v = get_dep_version([PY3, spades, '--version'], r'^.*SPAdes.*v([\d\.]*)')
    Log.msg('SPAdes is available:', v + ' ' + spades)

    return spades


# Bowtie 2
def dep_check_bowtie2(dir_dep, os_id, force):
    if os_id == 'mac':
        url = ('https://sourceforge.net/projects/bowtie-bio/files/bowtie2/'
               '2.4.5/bowtie2-2.4.5-macos-x86_64.zip/download')
    elif os_id == 'linux':
        url = ('https://sourceforge.net/projects/bowtie-bio/files/bowtie2/'
               '2.4.5/bowtie2-2.4.5-linux-x86_64.zip/download')

    dnld_path = opj(dir_dep, 'bowtie2.zip')

    try:
        if force is True:
            raise
        bowtie2 = which('bowtie2')
        bowtie2_build = which('bowtie2-build')
        run([bowtie2, '-h'])
        run([bowtie2_build, '-h'])
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'bowtie2'))
            bowtie2 = opj(dir_bin, 'bowtie2')
            bowtie2_build = opj(dir_bin, 'bowtie2-build')
            run([bowtie2, '-h'])
            run([bowtie2_build, '-h'])
        except Exception:
            Log.wrn('Bowtie 2 was not found on this system, trying to '
                    'download.')
            download_file(url, dnld_path)
            zip_ref = zipfile.ZipFile(dnld_path, 'r')
            zip_ref.extractall(dir_dep)
            zip_ref.close()

            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'bowtie2'))
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
                chmod(bowtie2 + bt2exe, stat.S_IRWXU | stat.S_IRGRP |
                      stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)

            if not ope(bowtie2):
                Log.err('Could not download Bowtie 2.')
                return None, None

    regexp = r'^.*?version\s([\d\.]*)'
    v = get_dep_version([bowtie2, '--version'], regexp)
    Log.msg('bowtie2 is available:', v + ' ' + bowtie2)
    v = get_dep_version([bowtie2_build, '--version'], regexp)
    Log.msg('bowtie2-build is available:', v + ' ' + bowtie2_build)

    return bowtie2, bowtie2_build


# Rcorrector
def dep_check_rcorrector(dir_dep, force):
    url = 'https://github.com/karolisr/Rcorrector/archive/master.tar.gz'
    dnld_path = opj(dir_dep, 'rcorrector.tar.gz')

    try:
        try:
            jellyfish = which('jellyfish')
            run([jellyfish, '--help'])
        except Exception:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'Rcorrector'))
            jellyfish = opj(dir_bin, 'jellyfish', 'bin', 'jellyfish')
            raise
        if force is True:
            raise
        rcorrector = which('run_rcorrector.pl')
        run([rcorrector, '-version'])
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'Rcorrector'))
            try:
                rcorrector = opj(dir_bin, 'run_rcorrector.pl')
                run([rcorrector, '-version'])
            except Exception:
                Log.wrn('Rcorrector was not found on this system, trying to '
                        'download.')
                raise
            try:
                run([jellyfish, '--version'])
            except Exception:
                Log.wrn(
                    'jellyfish is required by Rcorrector, but was not found. '
                    'Trying to download and recompile Rcorrector and '
                    'jellyfish.')
                raise
        except Exception:
            if ope(dnld_path):
                remove(dnld_path)
            if dir_bin != opj(dir_dep, ''):
                rmtree(dir_bin)
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(dir_dep)
            tar_ref.close()
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'Rcorrector'))
            try:
                Log.wrn('Compiling Rcorrector.')
                run('make', cwd=dir_bin)
                rcorrector = opj(dir_bin, 'run_rcorrector.pl')
                jellyfish = opj(dir_bin, 'jellyfish', 'bin', 'jellyfish')
                chmod(rcorrector, stat.S_IRWXU | stat.S_IRGRP |
                      stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
                run([rcorrector, '-version'])
                if not ope(jellyfish):
                    jellyfish = which('jellyfish')
                run([jellyfish, '--version'])
            except Exception:
                Log.err('Something went wrong while trying to compile '
                        'Rcorrector.')
                Log.msg('Try downloading and installing it manually from: '
                        'https://github.com/karolisr/Rcorrector')
                return None

    v = get_dep_version([rcorrector, '-version'], r'^Rcorrector\sv([\d\.]*)')
    Log.msg('Rcorrector is available:', v + ' ' + rcorrector)

    return rcorrector


# Kraken2
def dep_check_kraken2(dir_dep, os_id, release_name, force):
    url = 'https://github.com/karolisr/kraken2/archive/master.tar.gz'

    dnld_path = opj(dir_dep, 'kraken2.tar.gz')

    try:
        if force is True:
            raise
        kraken2 = which('kraken2')
        kraken2_build = which('kraken2-build')

        dir_bin = dirname(kraken2)
        classify_bin = opj(dir_bin, 'classify')
        _ = run([classify_bin], do_not_raise=True)
        if not _.stderr.startswith('classify: mandatory filename'):
            raise

        run([kraken2, '--help'])
        run([kraken2_build, '--help'])
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'kraken2'))
            kraken2 = opj(dir_bin, 'bin', 'kraken2')
            kraken2_build = opj(dir_bin, 'bin', 'kraken2-build')

            classify_bin = opj(dir_bin, 'bin', 'classify')
            _ = run([classify_bin], do_not_raise=True)
            if not _.stderr.startswith('classify: mandatory filename'):
                raise

            run([kraken2, '--help'])
            run([kraken2_build, '--help'])
        except Exception:
            Log.wrn('Kraken2 was not found on this system, trying to '
                    'download.')

            if ope(dnld_path):
                remove(dnld_path)

            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(dir_dep)
            tar_ref.close()

            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'kraken2'))
            classify_bin = opj(dir_bin, 'bin', 'classify')
            kraken2 = opj(dir_bin, 'bin', 'kraken2')
            kraken2_build = opj(dir_bin, 'bin', 'kraken2-build')

            makefile = opj(dir_bin, 'src', 'Makefile')
            replace_line_in_file(makefile,
                                 'cp $(PROGS) $(KRAKEN2_DIR)/',
                                 'cp $(PROGS) "$(KRAKEN2_DIR)"/')
            try:
                Log.wrn('Compiling Kraken2 Attempt 1')
                run(['./install_kraken2.sh', 'bin'], cwd=dir_bin)

                _ = run([classify_bin], do_not_raise=True)
                if not _.stderr.startswith('classify: mandatory filename'):
                    raise

                run([kraken2, '--help'])
                run([kraken2_build, '--help'])

            except Exception:
                try:
                    Log.wrn('Compiling Kraken2 Attempt 2')

                    dir_libomp = opj(dir_dep, 'libomp')

                    if ope(dir_libomp):
                        rmtree(dir_libomp)

                    libomp_fp, v = brew_get('libomp', os_id, release_name,
                                            dir_dep)

                    tar_ref = tarfile.open(libomp_fp, 'r:gz')
                    tar_ref.extractall(dir_dep)
                    tar_ref.close()

                    dir_libomp_l = opj(dir_libomp, v, 'lib')
                    dir_libomp_i = opj(dir_libomp, v, 'include')

                    if os_id == 'mac':
                        # Changes the shared library identification name of a
                        # dynamic shared library.
                        dylib_f = opj(dir_libomp_l, 'libomp.dylib')

                        chmod(dylib_f, stat.S_IRWXU | stat.S_IRUSR |
                              stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP |
                              stat.S_IROTH | stat.S_IWOTH)

                        cmd = ['install_name_tool', '-id', dylib_f, dylib_f]
                        run(cmd)

                        cxx_flags = ('CXXFLAGS = -L{} -I{} -Xpreprocessor '
                                     '-fopenmp -lomp -Wall -std=c++11 -O3')

                    elif os_id == 'linux':
                        cxx_flags = ('CXXFLAGS = -L{} -I{} -fopenmp -lomp '
                                     '-static -Wall -std=c++11 -O3')

                    cxx_flags = cxx_flags.format(dir_libomp_l, dir_libomp_i)

                    makefile = opj(dir_bin, 'src', 'Makefile')

                    replace_line_in_file(makefile,
                                         'CXXFLAGS = -fopenmp -Wall -std=c++11'
                                         ' -O3', cxx_flags)

                    run(['./install_kraken2.sh', 'bin'], cwd=dir_bin)

                    _ = run([classify_bin], do_not_raise=True)
                    if not _.stderr.startswith('classify: mandatory filename'):
                        raise

                    run([kraken2, '--help'])
                    run([kraken2_build, '--help'])

                except Exception:
                    try:
                        Log.wrn('Compiling Kraken2 Attempt 3')
                        makefile = opj(dir_bin, 'src', 'Makefile')
                        replace_line_in_file(makefile, cxx_flags,
                                             'CXXFLAGS = -Wall -std=c++11 -O3')
                        run(['./install_kraken2.sh', 'bin'], cwd=dir_bin)

                        _ = run([classify_bin], do_not_raise=True)
                        if not _.stderr.startswith('classify: mandatory filename'):
                            raise

                        run([kraken2, '--help'])
                        run([kraken2_build, '--help'])
                    except Exception:
                        pass

            if not ope(kraken2):
                Log.err('Something went wrong while trying to compile '
                        'Kraken2.')
                Log.msg('Try downloading and installing it manually from: '
                        'https://github.com/karolisr/kraken2')
                return None, None

    regexp = r'^.*?version\s([\d\.\-A-Za-z]*)'
    v = get_dep_version([kraken2, '--version'], regexp)
    Log.msg('kraken2 is available:', v + ' ' + kraken2)
    v = get_dep_version([kraken2_build, '--version'], regexp)
    Log.msg('kraken2-build is available:', v + ' ' + kraken2_build)

    return kraken2, kraken2_build


# ToDo: Refactor download_kraken2_dbs

def download_kraken2_dbs(dbs_path):
    base_kraken2_url = 'ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/'
    msg_prefix = 'Downloading Kraken2 database: '

    # ------------------------------------------------------------------------

    base = '16S_Silva132_20200326'
    url = base_kraken2_url + base + '.tgz'
    tgz = opj(dbs_path, base + '.tgz')
    # ToDo: Use pattern matching for the directory name.
    #       Instead of using 16S_SILVA132_k2db -> 16S_SILVA132
    p_orig = opj(dbs_path, '16S_SILVA132_k2db')
    db_name = '16S_Silva132'
    p_new = opj(dbs_path, db_name)

    if not ope(p_new):
        Log.msg(msg_prefix + db_name)
        download_file(url=url, local_path=tgz, protocol='ftp')
        tar_ref = tarfile.open(tgz, 'r:gz')
        tar_ref.extractall(dbs_path)
        tar_ref.close()
        remove(tgz)
        move(p_orig, p_new)

    # ------------------------------------------------------------------------

    base = '16S_Silva138_20200326'
    url = base_kraken2_url + base + '.tgz'
    tgz = opj(dbs_path, base + '.tgz')
    # ToDo: Use pattern matching for the directory name.
    #       Instead of using 16S_SILVA138_k2db -> 16S_SILVA138
    p_orig = opj(dbs_path, '16S_SILVA138_k2db')
    db_name = '16S_Silva138'
    p_new = opj(dbs_path, db_name)

    if not ope(p_new):
        Log.msg(msg_prefix + db_name)
        download_file(url=url, local_path=tgz, protocol='ftp')
        tar_ref = tarfile.open(tgz, 'r:gz')
        tar_ref.extractall(dbs_path)
        tar_ref.close()
        remove(tgz)
        move(p_orig, p_new)

    # ------------------------------------------------------------------------

    base = 'minikraken_8GB_202003'
    url = base_kraken2_url + base + '.tgz'
    tgz = opj(dbs_path, base + '.tgz')
    # ToDo: Use pattern matching for the directory name.
    #       Instead of using minikraken_8GB_20200312 -> minikraken
    p_orig = opj(dbs_path, 'minikraken_8GB_20200312')
    db_name = 'minikraken_8GB_2020-03-12'
    p_new = opj(dbs_path, db_name)

    if not ope(p_new):
        Log.msg(msg_prefix + db_name)
        download_file(url=url, local_path=tgz, protocol='ftp')
        tar_ref = tarfile.open(tgz, 'r:gz')
        tar_ref.extractall(dbs_path)
        tar_ref.close()
        remove(tgz)
        move(p_orig, p_new)

    # ------------------------------------------------------------------------

    base_dropbox_url = 'https://www.dropbox.com/s/'

    base = 'mitochondrion_and_plastid'
    garb = 'vkbp7iys6s76tvf/'
    url = base_dropbox_url + garb + base + '.tar.gz?dl=1'
    tgz = opj(dbs_path, base + '.tar.gz')
    p = opj(dbs_path, base)

    if not ope(p):
        Log.msg(msg_prefix + base)
        download_file(url=url, local_path=tgz, protocol='http')
        tar_ref = tarfile.open(tgz, 'r:gz')
        tar_ref.extractall(dbs_path)
        tar_ref.close()
        remove(tgz)

    # ------------------------------------------------------------------------

    base_dropbox_url = 'https://www.dropbox.com/s/'

    base = 'mitochondrion'
    garb = '6liwneb26uvjuec/'
    url = base_dropbox_url + garb + base + '.tar.gz?dl=1'
    tgz = opj(dbs_path, base + '.tar.gz')
    p = opj(dbs_path, base)

    if not ope(p):
        Log.msg(msg_prefix + base)
        download_file(url=url, local_path=tgz, protocol='http')
        tar_ref = tarfile.open(tgz, 'r:gz')
        tar_ref.extractall(dbs_path)
        tar_ref.close()
        remove(tgz)

    # ------------------------------------------------------------------------

    base_dropbox_url = 'https://www.dropbox.com/s/'

    base = 'plastid'
    garb = 's9vdg4mxrfy1szn/'
    url = base_dropbox_url + garb + base + '.tar.gz?dl=1'
    tgz = opj(dbs_path, base + '.tar.gz')
    p = opj(dbs_path, base)

    if not ope(p):
        Log.msg(msg_prefix + base)
        download_file(url=url, local_path=tgz, protocol='http')
        tar_ref = tarfile.open(tgz, 'r:gz')
        tar_ref.extractall(dbs_path)
        tar_ref.close()
        remove(tgz)

    # ------------------------------------------------------------------------

    base_dropbox_url = 'https://www.dropbox.com/s/'

    base = 'viral'
    garb = '7xz31c7vw088n27/'
    url = base_dropbox_url + garb + base + '.tar.gz?dl=1'
    tgz = opj(dbs_path, base + '.tar.gz')
    p = opj(dbs_path, base)

    if not ope(p):
        Log.msg(msg_prefix + base)
        download_file(url=url, local_path=tgz, protocol='http')
        tar_ref = tarfile.open(tgz, 'r:gz')
        tar_ref.extractall(dbs_path)
        tar_ref.close()
        remove(tgz)

    # ------------------------------------------------------------------------

    dbs_available, e = list_of_dirs_at_path(dbs_path)
    kraken2_dbs = {basename(p): p for p in dbs_available}
    return kraken2_dbs
