"""Differences between operating systems."""

import sys
from platform import mac_ver

import distro


def check_os():
    """Determine current operating system and distribution."""
    dist_id = None
    release_id = None
    release_name = ''

    os_id = None
    os_str = None

    if sys.platform == 'darwin':
        os_id = 'mac'
        os_str = 'macOS'
        release_id = mac_ver()[0]
        mv = release_id.split('.')
        if len(mv) == 3:
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
                elif mv[1] == '1':
                    release_name = 'Yosemite'
                else:
                    release_name = ''
            elif mv[0] == '11':
                if mv[1] == '0':
                    release_name = 'Big Sur'
                else:
                    release_name = ''

    elif sys.platform == 'win32':
        os_id = 'windows'
        os_str = 'Windows'

    elif sys.platform.startswith('linux'):
        os_id = 'linux'
        dist_id = distro.id()
        dist_name = distro.distro_release_attr('name')
        if dist_id == 'ubuntu' and dist_name == '':
            dist_name = dist_id.capitalize()
        os_str = 'Linux ({dist_name})'.format(dist_name=dist_name)

    debian_dists = ['debian', 'ubuntu']
    redhat_dists = ['centos', 'fedora', 'rhel', 'scientific']
    supported_dists = debian_dists + redhat_dists

    return {'os_id': os_id,
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
