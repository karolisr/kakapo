# -*- coding: utf-8 -*-

"""
Accounts for differences between operating systems.
"""

import distro
import sys

DEBIAN_DISTS = ['debian', 'ubuntu']
REDHAT_DISTS = ['centos', 'fedora', 'rhel', 'scientific']


def check_os():
    """
    Determine current operating system and distribution.
    """
    dist_id = None

    if sys.platform == 'darwin':
        os_id = 'mac'
        os_str = 'macOS'

    elif sys.platform == 'win32':
        os_id = 'windows'
        os_str = 'Windows'

        # print('Windows is not supported yet.')
        sys.exit(1)

    elif sys.platform.startswith('linux'):
        os_id = 'linux'
        dist_id = distro.id()
        dist_name = distro.distro_release_attr('name')
        if dist_id == 'ubuntu' and dist_name == '':
            dist_name = dist_id.capitalize()
        os_str = 'Linux ({dist_name})'.format(dist_name=dist_name)

        supported_dists = DEBIAN_DISTS + REDHAT_DISTS

        if dist_id not in supported_dists:
            # print('{dist_name} is not supported yet.'.format(
            #     dist_name=dist_name))
            sys.exit(1)

    else:
        # print('{p} is not supported yet.'.format(p=sys.platform))
        sys.exit(1)

    return os_id, os_str, dist_id

#
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
