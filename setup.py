"""Kakapo setup."""

# This is a workaround. pip, in some cases, will not find dependencies
# on GitHub, so we install them directly.
import sys
import subprocess

git_ntaxlocal = 'git+https://github.com/karolisr/ncbi-taxonomy-local'

if '--user' in sys.argv:
    subprocess.run([sys.executable, '-m', 'pip', 'install', '--upgrade',
                    '--user', git_ntaxlocal], check=False)
else:
    subprocess.run([sys.executable, '-m', 'pip', 'install', '--upgrade',
                    git_ntaxlocal], check=False)
# end pip workaround.

import setuptools  # noqa

from kakapo import __author__ as kakapo_author
from kakapo import __author_email__ as kakapo_author_email
from kakapo import __description__ as kakapo_description
from kakapo import __script_name__ as kakapo_script_name
from kakapo import __version__ as kakapo_version

with open('requirements.txt', 'r') as f:
    reqs = f.read()

reqs = reqs.split('\n')[0:-2]

setuptools.setup(
    name=kakapo_script_name,
    version=kakapo_version,
    author=kakapo_author,
    author_email=kakapo_author_email,
    description=kakapo_description,
    long_description=kakapo_description,
    long_description_content_type='text/markdown',
    url='https://github.com/karolisr/kakapo',
    packages=setuptools.find_packages(),
    install_requires=reqs,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: cc-by-sa-4.0',
        'Operating System :: OS Independent',
    ],
    entry_points={'console_scripts': ['kakapo=kakapo.__main__:main']}
)
