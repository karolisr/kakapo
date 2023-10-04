"""Differences between operating systems."""

import os
import sys
from os import popen, sysconf
from platform import mac_ver, machine

import distro


def check_os():
    """Determine current operating system and distribution."""
    dist_id = None
    release_id = None
    release_name = ''
    machine_type = machine()

    os_id = None
    os_str = None

    if sys.platform == 'darwin':
        os_id = 'mac'
        os_str = 'macOS'
        release_id = mac_ver()[0]
        mv = release_id.split('.')
        if len(mv) >= 2:
            if mv[0] == '10':
                if mv[1] == '15':
                    release_name = 'Catalina'
                elif mv[1] == '14':
                    release_name = 'Mojave'
                elif mv[1] == '13':
                    release_name = 'High Sierra'
                elif mv[1] == '12':
                    release_name = 'Sierra'
                elif mv[1] == '11':
                    release_name = 'El Capitan'
                elif mv[1] == '10':
                    release_name = 'Yosemite'
                else:
                    release_name = ''
            elif mv[0] == '11':
                release_name = 'Big Sur'
            elif mv[0] == '12':
                release_name = 'Monterey'
            elif mv[0] == '13':
                release_name = 'Ventura'
            elif mv[0] == '14':
                release_name = 'Sonoma'

    elif sys.platform == 'win32':
        os_id = 'windows'
        os_str = 'Windows'

    elif sys.platform.startswith('linux'):
        os_id = 'linux'
        dist_id = distro.id()
        dist_name = distro.distro_release_attr('name')
        dist = distro.linux_distribution()
        dist_name = dist[0]
        dist_version = dist[1]
        release_id = dist_version
        release_name = os_id
        if dist_id == 'ubuntu' and dist_name == '':
            dist_name = dist_id.capitalize()
        os_str = '{dist_name} {dist_version}'.format(dist_name=dist_name,
                                                     dist_version=dist_version)

    debian_dists = ['debian', 'ubuntu']
    redhat_dists = ['centos', 'fedora', 'rhel', 'scientific']
    supported_dists = debian_dists + redhat_dists

    return {'machine_type': machine_type,
            'os_id': os_id,
            'os_str': os_str,
            'dist_id': dist_id,
            'release_id': release_id,
            'release_name': release_name,
            'debian_dists': debian_dists,
            'redhat_dists': redhat_dists,
            'supported_dists': supported_dists}

# https://github.com/nir0s/distro
# possible distro values
#
# ==============  =========================================
# Distro ID       Distribution
# ==============  =========================================
# "amazon"        Amazon Linux
# "arch"          Arch Linux
# "centos"        CentOS
# "cloudlinux"    CloudLinux OS
# "debian"        Debian
# "exherbo"       Exherbo Linux
# "fedora"        Fedora
# "freebsd"       FreeBSD
# "gentoo"        GenToo Linux
# "ibm_powerkvm"  IBM PowerKVM
# "kvmibm"        KVM for IBM z Systems
# "linuxmint"     Linux Mint
# "mageia"        Mageia
# "mandriva"      Mandriva Linux
# "netbsd"        NetBSD
# "openbsd"       OpenBSD
# "opensuse"      openSUSE
# "oracle"        Oracle Linux (and Oracle Enterprise Linux)
# "parallels"     Parallels
# "pidora"        Pidora
# "raspbian"      Raspbian
# "rhel"          RedHat Enterprise Linux
# "scientific"    Scientific Linux
# "slackware"     Slackware
# "sles"          SUSE Linux Enterprise Server
# "ubuntu"        Ubuntu
# "xenserver"     XenServer
# ==============  =========================================


def cpu_count():
    os_id = check_os()['os_id']
    n_cores = 0
    n_threads = 0
    if os_id == 'linux':
        import psutil
        n_cores = psutil.cpu_count(logical=False)
        n_threads = psutil.cpu_count(logical=True)
    elif os_id == 'mac':
        n_cores = os.cpu_count()
        n_threads = os.cpu_count()

    return (n_cores, n_threads)


def sys_ram():
    os_id = check_os()['os_id']
    try:
        page_size = sysconf('SC_PAGE_SIZE')
        page_count = sysconf('SC_PHYS_PAGES')
        if page_size < 0 or page_count < 0:
            raise SystemError
        ram_b = sysconf('SC_PAGE_SIZE') * sysconf('SC_PHYS_PAGES')

    except ValueError:
        if os_id == 'mac':
            ram_b = int(
                float(popen("sysctl hw.memsize").readlines()[0].split()[1]))
        elif os_id == 'linux':
            ram_b = int(float(popen("free").readlines()[1].split()[1]) * 1024)
        else:
            raise NotImplementedError

    ram_g = ram_b / (1024 ** 3)
    return ram_g
