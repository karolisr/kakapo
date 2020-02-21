# kakapo [![Build Status](https://travis-ci.com/karolisr/kakapo.svg?branch=master)](https://travis-ci.com/karolisr/kakapo) [![License: cc-by-sa-4.0](https://i.creativecommons.org/l/by-sa/4.0/80x15.png)](http://creativecommons.org/licenses/by-sa/4.0/)

## Brief Description

Studies in many fields within life sciences increasingly rely on RNA sequencing (RNA-seq) data. As a result, RNA-seq datasets deposited to the NCBI Sequence Read Archive (SRA) are proliferating. In addition to serving as an archive for the original studies, these datasets present an opportunity for novel research.

Kakapo is a pipeline that allows users to extract and assemble a specified gene or a protein family from any number of SRA accessions (or their own RNA-seq data). Kakapo identifies open reading frames in the assembled transcripts and annotates them using InterProScan. Additionally, raw reads can be filtered for ribosomal, plastid, and mitochondrial reads or reads belonging to non-target organisms (viral, bacterial, etc.)

Kakapo can be flexibly employed to extract arbitrary loci, such as those commonly used for phylogenetic inference in systematics.

## Installation

Kakapo was designed for, and should work on, machines running macOS or Linux. It should be possible to use Kakapo on the [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) but I have not tested it. If you choose to try running Kakapo on Windows Subsystem for Linux, I suggest using the latest available Ubuntu distribution.

Kakapo supports Python 3 and will **not** work with Python 2. Use `pip` command below to install. In case you have both Python 2 and Python 3 on your system, or if you are not certain, make sure you have a Python 3 version of `pip` by running the command below:

```bash
pip -V
```

This should print the version of `pip` you have and if `pip` is using Python 3. If the output lists Python version 3.6 or higher, you are set.

```
pip 20.0.2 from .../python3.8/site-packages/pip (python 3.8)
```

Otherwise you may try:

```bash
pip3 -V
```

If `pip3` command works but `pip` does not, replace `pip` with `pip3` in the commands below:

```bash
pip install --upgrade git+https://github.com/karolisr/kakapo
```

In case this fails you may try:

```bash
pip install --user --upgrade git+https://github.com/karolisr/kakapo
```

or:

```bash
sudo -H pip install --upgrade git+https://github.com/karolisr/kakapo
```

In case none of the above commands work, you may not have `pip` installed. One easy way to install `pip` is by installing [Conda or Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).
