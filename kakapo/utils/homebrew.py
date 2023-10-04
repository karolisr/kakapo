"""Homebrew API."""

from os.path import abspath, expanduser
from os.path import join as opj

from kakapo.utils.http import download_file, get


def brew_get(package, os, platform, machine_type, dnld_dir):

    url_mac = 'https://formulae.brew.sh/api/formula/{}.json'
    url_linux = 'https://formulae.brew.sh/api/formula-linux/{}.json'

    if os == 'mac':
        if platform == 'Ventura':
            if machine_type == 'x86_64':
                platform = 'ventura'
            else:
                platform = 'arm64_ventura'
        elif platform == 'Monterey':
            if machine_type == 'x86_64':
                platform = 'monterey'
            else:
                platform = 'arm64_monterey'
        elif platform == 'Big Sur':
            if machine_type == 'x86_64':
                platform = 'big_sur'
            else:
                platform = 'arm64_big_sur'
        elif platform == 'Catalina':
            platform = 'catalina'
        elif platform == 'Mojave':
            platform = 'mojave'
        elif platform == 'High Sierra':
            platform = 'high_sierra'
        elif platform == 'Sierra':
            platform = 'sierra'
        elif platform == 'El Capitan':
            platform = 'el_capitan'
        elif platform == 'Yosemite':
            platform = 'yosemite'
        else:
            return None, None

        url = url_mac.format(package)

    elif os == 'linux':
        if platform == 'linux':
            platform = 'x86_64_linux'
        else:
            return None, None

        url = url_linux.format(package)

    else:
        return None, None

    try:
        r = get(url)
    except Exception as e:
        print(e)
        return None, None

    data = r.json()

    version = data['versions']['stable']

    dnld_dir = abspath(expanduser(dnld_dir))
    file_url = data['bottle']['stable']['files'][platform]['url']
    # file_name = file_url.split('/')[-1]
    file_name = '{}__{}__{}__{}.tar.gz'.format(package, version, os, platform)
    dnld_path = opj(dnld_dir, file_name)

    # https://github.com/Homebrew/brew/pull/11070#issuecomment-830448946
    download_file(file_url, dnld_path, {'Authorization': 'Bearer QQ=='})

    return dnld_path, version
