"""Homebrew API."""

from os.path import abspath
from os.path import expanduser
from os.path import join as opj

from kakapo.utils.http import get
from kakapo.utils.http import download_file


def brew_get(package, os, platform, dnld_dir):

    url_mac = 'https://formulae.brew.sh/api/formula/{}.json'
    url_linux = 'https://formulae.brew.sh/api/formula-linux/{}.json'

    if os == 'mac':
        if platform == 'Catalina':
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
            return None

        url = url_mac.format(package)

    if os == 'linux':
        if platform == 'linux':
            platform = 'x86_64_linux'
        else:
            return None

        url = url_linux.format(package)

    r = get(url)

    try:
        r.raise_for_status()
    except Exception as e:
        print(e)
        return None

    data = r.json()

    version = data['versions']['stable']

    dnld_dir = abspath(expanduser(dnld_dir))
    file_url = data['bottle']['stable']['files'][platform]['url']
    file_name = file_url.split('/')[-1]
    dnld_path = opj(dnld_dir, file_name)

    download_file(file_url, dnld_path)

    return dnld_path, version
