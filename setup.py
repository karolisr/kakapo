import setuptools  # noqa

with open('README.md', 'r') as f:
    long_description = f.read()

description = 'Extract and annotate protein family members from transcriptomes.'

with open('requirements.txt', 'r') as f:
    reqs = f.read()

reqs = reqs.split('\n')[0:-2]

setuptools.setup(
    name='kakapo',
    version='0.0.1',
    author='Karolis Ramanauskas',
    author_email='kraman2@uic.edu',
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/karolisr/kakapo',
    packages=setuptools.find_packages(),
    install_requires=reqs + ['ncbi-taxonomy-local @ https://github.com/karolisr/ncbi-taxonomy-local/archive/master.zip'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: cc-by-sa-4.0',
        'Operating System :: OS Independent',
    ],
    entry_points={'console_scripts': ['kakapo=kakapo.__main__:main']}
)
