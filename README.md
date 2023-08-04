# kakapo [![Build Status](https://app.travis-ci.com/karolisr/kakapo.svg?branch=master)](https://app.travis-ci.com/karolisr/kakapo)

## Brief Description

Studies in many fields within life sciences increasingly rely on RNA sequencing (RNA-seq) data. As a result, RNA-seq datasets deposited to the NCBI Sequence Read Archive (SRA) are proliferating. In addition to serving as an archive for the original studies, these datasets present an opportunity for novel research.

Kakapo is a pipeline that allows users to extract and assemble a specified gene or a protein family from any number of SRA accessions (or their own RNA-seq data). Kakapo identifies open reading frames in the assembled transcripts and annotates them using InterProScan. Additionally, raw reads can be filtered for ribosomal, plastid, and mitochondrial reads or reads belonging to non-target organisms (viral, bacterial, etc.)

Kakapo can be flexibly employed to extract arbitrary loci, such as those commonly used for phylogenetic inference in systematics.

Brief overview of kakapo from the Botany 2020 conference: [https://youtu.be/2D04DQlV6CA](https://youtu.be/2D04DQlV6CA)

## Installation

See [Wiki](https://github.com/karolisr/kakapo/wiki) for the [installation instructions](https://github.com/karolisr/kakapo/wiki/Installation).
