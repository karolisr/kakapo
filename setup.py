"""Kakapo Setup Script."""

from setuptools import setup
from setuptools import find_packages

from kakapo import __author__ as kakapo_author
from kakapo import __author_email__ as kakapo_author_email
from kakapo import __description__ as kakapo_description
from kakapo import __script_name__ as kakapo_script_name
from kakapo import __url__ as kakapo_url
from kakapo import __version__ as kakapo_version

with open('requirements.txt', 'r') as f:
    _ = f.read()

reqs = _.splitlines()

setup(name=kakapo_script_name,
      version=kakapo_version,
      description=kakapo_description,
      long_description=kakapo_description,
      long_description_content_type='text/markdown',
      author=kakapo_author,
      author_email=kakapo_author_email,
      maintainer=kakapo_author,
      maintainer_email=kakapo_author_email,
      url=kakapo_url,
      include_package_data=True,
      install_requires=reqs,
      packages=find_packages(),
      python_requires='>=3.9',
      entry_points={'console_scripts': ['kakapo=kakapo.__main__:run_kakapo']},
      classifiers=[
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'Natural Language :: English',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: POSIX :: Linux',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3.11',
          'Programming Language :: Python :: 3.12',
          'Programming Language :: Python :: 3.13',
          'Programming Language :: Python :: 3 :: Only',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ]
      )
