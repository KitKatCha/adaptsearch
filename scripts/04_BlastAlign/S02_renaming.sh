#!/bin/bash
# Author : Victor Mataigne (ABiMS, Station Biologique de Roscoff)
# A small script which rename a file correctly, according to its number of sequences

name=$(echo $1 | cut -d _ -f 1)
under=$(echo _sp)
number=$(grep '/' $1 | wc -l)  
phy=$(echo .phy)
new_file=$name$under$number$phy
cat $1 > $new_file
#mv $1 $new_file
#ln -s $new_file $1

name=$(echo $2 | cut -d _ -f 1)
under=$(echo _sp)   
number=$(head -4 $2 | tail -1 | cut -d " "  -f 2 | cut -d = -f 2)    
nxs=$(echo .nxs)
new_file=$name$under$number$nxs
cat $2 > new_file  
#mv $2 $new_file
#ln -s $new_file $2

if [ ! -z "$3" ]
then
    name=$(echo $3 | cut -d _ -f 1)
    under=$(echo _sp)   
    number=$(grep '>' $3 | wc -l)
    fasta=$(echo .fasta)
    new_file=$name$under$number$fasta
    n=$(echo $3) 
    cat $3 > $new_file 
    #mv $3 $new_file
    #ln -s $new_file $n
fi
