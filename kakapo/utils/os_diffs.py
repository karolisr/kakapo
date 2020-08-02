"""Differences between operating systems."""

import sys

import distro


def check_os():
    """Determine current operating system and distribution."""
    dist_id = None

    os_id = None
    os_str = None

    if sys.platform == 'darwin':
        os_id = 'mac'
        os_str = 'macOS'

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
