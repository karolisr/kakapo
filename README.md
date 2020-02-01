# kakapo [![Build Status](https://travis-ci.com/karolisr/kakapo.svg?branch=master)](https://travis-ci.com/karolisr/kakapo) [![License: cc-by-sa-4.0](https://i.creativecommons.org/l/by-sa/4.0/80x15.png)](http://creativecommons.org/licenses/by-sa/4.0/)

Studies in many fields within life sciences increasingly rely on RNA sequencing (RNA-seq) data. As a result, RNA-seq datasets deposited to the NCBI Sequence Read Archive (SRA) are proliferating. In addition to serving as an archive for the original studies, these datasets present an opportunity for novel research.

Kakapo is a pipeline that allows users to extract and assemble a specified gene or a protein family from any number of SRA accessions (or their own RNA-seq data). Kakapo identifies open reading frames in the assembled transcripts and annotates them using InterProScan. Additionally, raw reads can be filtered for ribosomal, plastid, and mitochondrial reads or reads belonging to non-target organisms (viral, bacterial, etc.)

Kakapo can be flexibly employed to extract arbitrary loci, such as those commonly used for phylogenetic inference in systematics.

Installation:

```bash
pip install --upgrade git+https://github.com/karolisr/kakapo
```

In case this fails you may try:

```bash
sudo -H pip install --upgrade git+https://github.com/karolisr/kakapo
```

or:

```bash
pip install --user --upgrade git+https://github.com/karolisr/kakapo
```
