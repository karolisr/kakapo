import setuptools  # noqa

from kakapo import __version__ as kakapo_version
from kakapo import __script_name__ as kakapo_script_name
from kakapo import __author__ as kakapo_author
from kakapo import __author_email__ as kakapo_author_email
from kakapo import __description__ as kakapo_description

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
    install_requires=reqs + ['ncbi-taxonomy-local @ https://github.com/karolisr/ncbi-taxonomy-local/archive/master.zip'],  # noqa
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: cc-by-sa-4.0',
        'Operating System :: OS Independent',
    ],
    entry_points={'console_scripts': ['kakapo=kakapo.__main__:main']}
)
