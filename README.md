[![Build Status](https://travis-ci.org/Mataivic/adaptsearch.svg?branch=master)](https://travis-ci.org/Mataivic/adaptsearch)

# AdaptSearch: A Galaxy pipeline for sequence comparison and the search of positive selection from orthologous groups derived from RNAseq datasets

## Abstract:
This Galaxy workflow use as input transcriptome assemblies for n species and proceeds to find putative orthogroups. Multiple alignments are then produced with the corresponding full-length transcripts within each orthogroup. After finding the coding sequences in the correct frame and indels removal, the suite reconstructs a phylogenomic tree from a concatenated set of the coding-sequence alignments, and search for positively-selected genes along the branches of the tree and positively-selected codons along both the whole sequence alignment or each transcript, separately. In addition, AdaptSearch estimates base, codon and amino-acid frequencies from the selected set of genes and their products as well as all the substitu-tions from one species to another and tests a priori ecological hypotheses using a codon resampling pro-cedure or binomial tests on the sequence composition and/or dN and dS substitution rates.

## Availability:
### Ready to use
AdaptSearch can be used on [galaxy.sb-roscoff.fr](http://galaxy.sb-roscoff.fr/). An account is require to access to the Galaxy instance: http://abims.sb-roscoff.fr/account

### Installation on a Galaxy instance
The tools can be installed in any Galaxy instance using the Galaxy Tool Shed: [abims-sbr/suite_adaptsearch/](https://toolshed.g2.bx.psu.edu/view/abims-sbr/suite_adaptsearch/)

## Example Dataset
A toy example to test the pipeline is available on Zenodo: https://zenodo.org/record/3356025#.XUFcDJMzZhE 

## Other repositories
[Code and user documentation](https://github.com/Mataivic/AdaptSearch_Documentation)

[Visualization tools](https://github.com/Mataivic/AdaptSearch_visualization_tools)

## Credits
- Victor Mataigne
[ABiMS - Station Biologique de Roscoff - France - CNRS/SU](http://abims.sb-roscoff.fr/)
[ABICE / Adaptation et Biologie des Invertébrés en Conditions Extrêmes](http://www.sb-roscoff.fr/fr/abice-adaptation-et-biologie-des-invertebres-en-conditions-extremes)
- Eric Fontanillas
[ABICE / Adaptation et Biologie des Invertébrés en Conditions Extrêmes](http://www.sb-roscoff.fr/fr/abice-adaptation-et-biologie-des-invertebres-en-conditions-extremes)
- Julie Baffard
[ABiMS - Station Biologique de Roscoff - France - CNRS/SU](http://abims.sb-roscoff.fr/)
- Misharl Monsoor
[ABiMS - Station Biologique de Roscoff - France - CNRS/SU](http://abims.sb-roscoff.fr/)
- Gildas Le Corguillé
[ABiMS - Station Biologique de Roscoff - France - CNRS/SU](http://abims.sb-roscoff.fr/)
- Didier Jollivet
[ABICE / Adaptation et Biologie des Invertébrés en Conditions Extrêmes](http://www.sb-roscoff.fr/fr/abice-adaptation-et-biologie-des-invertebres-en-conditions-extremes)


