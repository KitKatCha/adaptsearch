#!/usr/bin/env python
# coding: utf8
## Author: Eric Fontanillas
## Modification: 03/09/14 by Julie BAFFARD
## Last modification : 05/03/18 by Victor Mataigne

## Description : find and remove indels

####################
###### DEF 2 #######
####################
def remove_position_with_too_much_missing_data(bash_aa, bash_nuc, MIN_SPECIES_NB):

    ## 1 ## Get alignment length
    fasta_name0 = bash_aa.keys()[0]
    ln_aa = len(bash_aa[fasta_name0])

    ln_nuc = len(bash_nuc[fasta_name0])


    ## 2 ## Get positions keeped in aa alignment
    LIST_POSITION_KEEPED_aa = []
    i=0
    while i < ln_aa:
        site = []
        for fasta_name in bash_aa.keys():
            pos = bash_aa[fasta_name][i]

            if pos != "-" and pos != "?" and pos != "X":
                site.append(pos)
        if len(site) >= MIN_SPECIES_NB:
            LIST_POSITION_KEEPED_aa.append(i)
        i = i+1

    ## 3 ## Get positions keeped in nuc alignment
    LIST_POSITION_KEEPED_nuc = []
    for position in LIST_POSITION_KEEPED_aa:
        position1 = position*3
        position2 = position*3 + 1
        position3 = position*3 + 2
        LIST_POSITION_KEEPED_nuc.append(position1)
        LIST_POSITION_KEEPED_nuc.append(position2)
        LIST_POSITION_KEEPED_nuc.append(position3)

    ## 4 ## Create entries for "filtered_bash" for aa & nuc
    filtered_bash_aa = {}
    filtered_bash_nuc = {}
    for fasta_name in bash_aa.keys():
        filtered_bash_aa[fasta_name] = ""
    for fasta_name in bash_nuc.keys():
        filtered_bash_nuc[fasta_name] = ""

    ## 5 ## Write "filtered_bash" for aa
    j=0
    while j < ln_aa:
        for fasta_name in bash_aa.keys():
            seq=filtered_bash_aa[fasta_name]
            pos=bash_aa[fasta_name][j]

            if j in LIST_POSITION_KEEPED_aa:
                seq = seq + pos
                filtered_bash_aa[fasta_name] = seq
        j = j + 1

    ## 6 ## Remove empty sequence
    for name in filtered_bash_aa.keys():
        seq = filtered_bash_aa[name]
        if seq == '':
            del filtered_bash_aa[name]


    ## 7 ## Write "filtered_bash" for nuc
    j=0
    while j < ln_nuc:
        for fasta_name in bash_nuc.keys():
            seq=filtered_bash_nuc[fasta_name]
            #print seq
            pos=bash_nuc[fasta_name][j]

            if j in LIST_POSITION_KEEPED_nuc:
                seq = seq + pos
                filtered_bash_nuc[fasta_name] = seq
        j = j + 1

    ## 8 ## Remove empty sequence
    for name in filtered_bash_nuc.keys():
        seq = filtered_bash_nuc[name]
        if seq == '':
            del filtered_bash_nuc[name]

    return(filtered_bash_aa, filtered_bash_nuc)
####################################


#######################
##### RUN RUN RUN #####
#######################
import string, os, time, re, sys
from dico import dico

### 0 ### PARAMETERS
MIN_SPECIES_NB = int(sys.argv[1])
MIN_LENGTH_FINAL_ALIGNMENT_NUC = int(sys.argv[2])
n0 = 0
bad = 0
good = 0
list_new_file = []
dicoco = {}
list_file = []
name_elems = ["orthogroup", "0", "with", "0", "species.fasta"]

### 1 ### IN
path_IN1 = "./07_CDS_aa/"
L_IN1 = os.listdir(path_IN1)
lenght = len(L_IN1)
path_IN2 = "./07_CDS_nuc/"
L_IN2 = os.listdir(path_IN2)

## 2 ## OUT
os.mkdir("08_CDS_aa_MINIMUM_MISSING_SEQUENCES")
path_OUT1 = "08_CDS_aa_MINIMUM_MISSING_SEQUENCES"
os.mkdir("08_CDS_nuc_MINIMUM_MISSING_SEQUENCES")
path_OUT2 = "08_CDS_nuc_MINIMUM_MISSING_SEQUENCES"


for file in L_IN1:
    file_INaa = "%s/%s" %(path_IN1, file)
    file_INnuc = "%s/%s" %(path_IN2, file)

    dico_aa = dico(file_INaa)   ### DEF 1 ###
    dico_nuc = dico(file_INnuc)   ### DEF 1 ###

    if len(dico_aa) < MIN_SPECIES_NB :
        list_file.append(file)

if list_file == lenght :
    MIN_SPECIES_NB == MIN_SPECIES_NB - 1


for file in L_IN1 :
    file_INaa = "%s/%s" %(path_IN1, file)
    file_INnuc = "%s/%s" %(path_IN2, file)

    dico_aa = dico(file_INaa)   ### DEF 1 ###
    dico_nuc = dico(file_INnuc)   ### DEF 1 ###

    ## 4.1 ## REMOVE POSITIONS WITH TOO MUCH MISSING DATA (i.e. not enough taxa represented at each position in the alignment)
    filtered_bash_aa, filtered_bash_nuc = remove_position_with_too_much_missing_data(dico_aa, dico_nuc, MIN_SPECIES_NB)   ### DEF 2 ###

    k = filtered_bash_nuc.keys()
    new_leng_nuc = 0
    if k != []:
        k0 = k[0]
        seq0 = filtered_bash_nuc[k0]
        new_leng_nuc = len(seq0)

    ## 4.3 ## Change file name for output, depending the number of species remaining in the alignment 
    n0+=1
    #name_elems[1] = str(n0)
    name_elems[1] = file.split('_')[1]
    name_elems[3] =  str(len(filtered_bash_aa.keys()))
    new_name = "_".join(name_elems)

    ## 4.5 ## Write filtered alignment in OUTPUTs
    ## aa
    if filtered_bash_aa != {} and new_leng_nuc >= MIN_LENGTH_FINAL_ALIGNMENT_NUC:
        OUTaa=open("%s/%s" %(path_OUT1, new_name), "w")
        for fasta_name in filtered_bash_aa.keys():
            seq_aa = filtered_bash_aa[fasta_name]
            OUTaa.write("%s\n" %fasta_name)
            OUTaa.write("%s\n" %seq_aa)
        OUTaa.close()
    # nuc
    if filtered_bash_nuc != {} and new_leng_nuc >= MIN_LENGTH_FINAL_ALIGNMENT_NUC:
        good+=1
        OUTnuc=open("%s/%s" %(path_OUT2, new_name), "w")
        for fasta_name in filtered_bash_nuc.keys():
            seq_nuc = filtered_bash_nuc[fasta_name]
            OUTnuc.write("%s\n" %fasta_name[0:3])
            OUTnuc.write("%s\n" %seq_nuc)
        OUTnuc.close()
    else:
        bad+=1


## 5 ## Print
print "*************** 2nd Filter : removal of the indel ***************"
print "\nTotal number of locus recorded  = %d" %n0
print "\tTotal number of locus with no indels (SAVED) = %d" %good
print "\tTotal number of locus, when removing indel, wich are empty (EXCLUDED) = %d" %bad
print ""
