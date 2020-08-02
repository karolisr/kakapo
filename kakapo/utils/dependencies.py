"""Install Kakapo dependencies."""

import stat
import tarfile
import zipfile
from os import chmod
from os import linesep as lns
from os import remove
from os.path import exists as ope
from os.path import join as opj
from shutil import move
from shutil import rmtree
from tempfile import NamedTemporaryFile

from kakapo.utils.logging import Log
from kakapo.utils.misc import download_file
from kakapo.utils.misc import list_of_dirs_at_path
from kakapo.utils.misc import replace_line_in_file
from kakapo.utils.subp import run
from kakapo.utils.subp import run_then_grep
from kakapo.utils.subp import which


def install_deps():
    Log.inf('Checking for dependencies.')

    seqtk = dep_check_seqtk()
    trimmomatic, adapters = dep_check_trimmomatic()
    fasterq_dump = dep_check_sra_toolkit()
    makeblastdb, _, tblastn = dep_check_blast()
    vsearch = dep_check_vsearch()
    spades = dep_check_spades()
    bowtie2, bowtie2_build = dep_check_bowtie2()
    rcorrector = dep_check_rcorrector()
    kraken2, kraken2_build = dep_check_kraken2()


def dnld_kraken2_dbs():
    Log.inf('Checking for available Kraken2 databases.')

    kraken2_dbs = download_kraken2_dbs()
    for db in sorted(kraken2_dbs.keys()):
        Log.msg('Found Kraken2 database:', db)


# ----------------------------------------------------------------------------


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
        run(cmd)
    except Exception:
        try:
            seqtk = opj(dir_bin, 'seqtk')
            cmd[0] = seqtk
            run(cmd)
        except Exception:
            Log.wrn('Seqtk was not found on this system, trying to download.')
            download_file(url, dnld_path)
            zip_ref = zipfile.ZipFile(dnld_path, 'r')
            zip_ref.extractall(dir_dep)
            zip_ref.close()
            try:
                Log.wrn('Compiling Seqtk.')
                run('make', cwd=dir_bin)
                run(cmd)
            except Exception:
                replace_line_in_file(opj(dir_bin, 'Makefile'),
                                     'CC=gcc', 'CC=cc')
                try:
                    run('make', cwd=dir_bin)
                    run(cmd)
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
                        r'^\s*(.+)\s*$')
    Log.msg('Trimmomatic is available:', v + ' ' + trimmomatic)

    path_adapters = _write_trimmomatic_adapters_file()

    return trimmomatic, path_adapters


# SRA Toolkit
def dep_check_sra_toolkit(dir_dep, os_id, dist_id, debian_dists, redhat_dists,
                          force):
    if os_id == 'mac':
        url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/'
               'sratoolkit.2.9.6-mac64.tar.gz')
    elif os_id == 'linux':
        if dist_id in debian_dists:
            url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.2/'
                   'sratoolkit.2.10.2-ubuntu64.tar.gz')
        elif dist_id in redhat_dists:
            url = ('https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.2/'
                   'sratoolkit.2.10.2-centos_linux64.tar.gz')

    dnld_path = opj(dir_dep, 'sra-toolkit.tar.gz')

    fasterq_dump = None
    try:
        if force is True:
            raise
        fasterq_dump = which('fasterq-dump')
        run(fasterq_dump)
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'sratoolkit'))
            fasterq_dump = opj(dir_bin, 'bin', 'fasterq-dump')
            run(fasterq_dump)
        except Exception:
            Log.wrn('SRA Toolkit was not found on this system, trying to '
                    'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(dir_dep)
            tar_ref.close()

            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'sratoolkit'))
            fasterq_dump = opj(dir_bin, 'bin', 'fasterq-dump')

            if not ope(fasterq_dump):
                Log.err('Could not download SRA Toolkit.')
                return None

    v = get_dep_version([fasterq_dump, '--version'], r':\s([\d\.]*)')
    Log.msg('fasterq-dump is available:', v + ' ' + fasterq_dump)

    return fasterq_dump


# BLAST+
def dep_check_blast(dir_dep, os_id, dist_id, debian_dists, redhat_dists,
                    force):
    if os_id == 'mac':
        url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/'
               'ncbi-blast-2.10.0+-x64-macosx.tar.gz')
    elif os_id == 'linux':
        if dist_id in debian_dists:
            url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/'
                   '2.10.0/ncbi-blast-2.10.0+-x64-linux.tar.gz')
        elif dist_id in redhat_dists:
            url = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/'
                   '2.10.0/ncbi-blast-2.10.0+-x64-linux.tar.gz')

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
def dep_check_vsearch(dir_dep, force):
    url = 'https://github.com/torognes/vsearch/archive/master.tar.gz'
    dnld_path = opj(dir_dep, 'vsearch.tar.gz')
    dir_bin = opj(dir_dep, 'vsearch-master')

    try:
        if force is True:
            raise
        vsearch = which('vsearch')
        run(vsearch)
    except Exception:
        try:
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
                Log.wrn('Compiling Vsearch.')
                run('./autogen.sh', cwd=dir_bin)
                run('./configure', cwd=dir_bin)
                run('make', cwd=dir_bin)
                vsearch = opj(dir_bin, 'bin', 'vsearch')
                run(vsearch)
            except Exception:
                Log.err(
                    'Something went wrong while trying to compile Vsearch.')
                Log.msg('Try downloading and installing it manually from: '
                        'https://github.com/torognes/vsearch')
                return None

    v = get_dep_version([vsearch, '-version'], r'vsearch\sv([\d\.]*)')
    Log.msg('Vsearch is available:', v + ' ' + vsearch)

    return vsearch


# SPAdes
def dep_check_spades(dir_dep, force, os_id):
    if os_id == 'mac':
        url = ('http://cab.spbu.ru/files/release3.14.0/'
               'SPAdes-3.14.0-Darwin.tar.gz')
    elif os_id == 'linux':
        url = ('http://cab.spbu.ru/files/release3.14.0/'
               'SPAdes-3.14.0-Linux.tar.gz')

    dnld_path = opj(dir_dep, 'SPAdes.tar.gz')

    try:
        if force is True:
            raise
        spades = which('spades.py')
        run(spades)
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'SPAdes'))
            spades = opj(dir_bin, 'bin', '../tools/spades.py')
            run(spades)
        except Exception:
            Log.wrn('SPAdes was not found on this system, trying to download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(dir_dep)
            tar_ref.close()
            try:
                dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'SPAdes'))
                spades = opj(dir_bin, 'bin', '../tools/spades.py')
                run(spades)
            except Exception:
                Log.err('Could not download SPAdes.')
                return None

    v = get_dep_version([spades, '--version'], r'^.*SPAdes.*v([\d\.]*)')
    Log.msg('SPAdes is available:', v + ' ' + spades)

    return spades


# Bowtie 2
def dep_check_bowtie2(dir_dep, force, os_id):
    if os_id == 'mac':
        url = ('https://sourceforge.net/projects/bowtie-bio/files/bowtie2/'
               '2.4.1/bowtie2-2.4.1-macos-x86_64.zip/download')
    elif os_id == 'linux':
        url = ('https://sourceforge.net/projects/bowtie-bio/files/bowtie2/'
               '2.4.1/bowtie2-2.4.1-linux-x86_64.zip/download')

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
                run([jellyfish, '--help'])
            except Exception:
                Log.wrn(
                    'jellyfish is required by Rcorrector, but was not found. '
                    'Trying to download and recompile Rcorrector and '
                    'jellyfish.')
                raise
        except Exception:
            if ope(dnld_path):
                remove(dnld_path)
            dir_no_rc = opj(dir_dep, '')
            if dir_bin != dir_no_rc:
                rmtree(dir_bin)
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(dir_dep)
            tar_ref.close()
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'Rcorrector'))
            # jellyfish ######################################################
            jf_tgz_path = opj(dir_bin, 'jellyfish.tar.gz')
            jf_url = ('https://github.com/gmarcais/Jellyfish/releases/'
                      'download/v2.3.0/jellyfish-2.3.0.tar.gz')
            download_file(jf_url, jf_tgz_path)
            tar_ref = tarfile.open(jf_tgz_path, 'r:gz')
            tar_ref.extractall(dir_bin)
            tar_ref.close()
            ##################################################################
            try:
                Log.wrn('Compiling Rcorrector.')
                run('make', cwd=dir_bin)
                rcorrector = opj(dir_bin, 'run_rcorrector.pl')
                chmod(rcorrector, stat.S_IRWXU | stat.S_IRGRP |
                      stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
                run([rcorrector, '-version'])
                run([jellyfish, '--help'])
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
def dep_check_kraken2(dir_dep, force):
    url = 'https://github.com/karolisr/kraken2/archive/master.tar.gz'

    dnld_path = opj(dir_dep, 'kraken2.tar.gz')

    try:
        if force is True:
            raise
        kraken2 = which('kraken2')
        kraken2_build = which('kraken2-build')
        run([kraken2, '--help'])
        run([kraken2_build, '--help'])
    except Exception:
        try:
            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'kraken2'))
            kraken2 = opj(dir_bin, 'bin', 'kraken2')
            kraken2_build = opj(dir_bin, 'bin', 'kraken2-build')
            run([kraken2, '--help'])
            run([kraken2_build, '--help'])
        except Exception:
            Log.wrn('Kraken2 was not found on this system, trying to '
                    'download.')
            download_file(url, dnld_path)
            tar_ref = tarfile.open(dnld_path, 'r:gz')
            tar_ref.extractall(dir_dep)
            tar_ref.close()

            dir_bin = opj(dir_dep, get_dep_dir(dir_dep, 'kraken2'))
            kraken2 = opj(dir_bin, 'bin', 'kraken2')
            kraken2_build = opj(dir_bin, 'bin', 'kraken2-build')

            makefile = opj(dir_bin, 'src', 'Makefile')
            replace_line_in_file(makefile,
                                 '\tcp $(PROGS) $(KRAKEN2_DIR)/',
                                 '\tcp $(PROGS) "$(KRAKEN2_DIR)"/')
            try:
                Log.wrn('Compiling Kraken2')
                run(['./install_kraken2.sh', 'bin'], cwd=dir_bin)
            except Exception:
                try:
                    makefile = opj(dir_bin, 'src', 'Makefile')
                    replace_line_in_file(makefile,
                                         'CXXFLAGS = -fopenmp -Wall -std=c++11'
                                         '-O3',
                                         'CXXFLAGS = -Wall -std=c++11 -O3')
                    run(['./install_kraken2.sh', 'bin'], cwd=dir_bin)
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

    base = '16S_Silva_20190418'
    url = base_kraken2_url + base + '.tgz'
    tgz = opj(dbs_path, base + '.tgz')
    p_orig = opj(dbs_path, base)
    db_name = '16S_Silva'
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

    base = 'minikraken2_v2_8GB_201904_UPDATE'
    url = base_kraken2_url + base + '.tgz'
    tgz = opj(dbs_path, base + '.tgz')
    p_orig = opj(dbs_path, base)
    db_name = 'minikraken2_v2'
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
    kraken2_dbs = {k: opj(dbs_path, k) for k in dbs_available}
    return kraken2_dbs
