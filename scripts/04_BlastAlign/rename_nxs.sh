#!/bin/bash
# Author : Victor Mataigne (ABiMS, Station Biologique de Roscoff)
# A small script which rename a file correctly, according to its number of sequences

name=$(echo $1 | cut -d _ -f 1)
under=$(echo _sp)   
number=$(head -4 $1 | tail -1 | cut -d " "  -f 2 | cut -d = -f 2)    
nxs=$(echo .nxs)
new_file=$name$under$number$nxs
cat $1 > $new_file
ln -s $new_file $2