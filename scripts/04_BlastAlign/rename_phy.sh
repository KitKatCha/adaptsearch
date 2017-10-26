#!/bin/bash
# Author : Victor Mataigne (ABiMS, Station Biologique de Roscoff)
# A small script which rename a file correctly, according to its number of sequences

# $1 : locusX_spY.fasta.phy
# $2 : locusX_spY.fasta.nxs

name=$(echo $1 | cut -d _ -f 1)
under=$(echo _sp)
number=$(grep '/' $1 | wc -l) 
phy=$(echo .phy)
new_file=$name$under$number".phy"
mv $1 new_file
new_file=$name$under$number".nxs"
mv $2 new_file