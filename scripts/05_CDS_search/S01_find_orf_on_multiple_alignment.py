#!/usr/bin/env python
# coding: utf8
# Author: Eric Fontanillas
# Modification: 03/09/14 by Julie BAFFARD
# Last modification : 25/07/18 by Victor Mataigne

# Description: Predict potential ORF on the basis of 2 criteria + 1 optional criteria
                # CRITERIA 1 - Longest part of the alignment of sequence without codon stop "*", tested in the 3 potential ORF
                # CRITERIA 2 - This longest part should be > 150nc or 50aa
                # CRITERIA 3 - [OPTIONNAL] A codon start "M" should be present in this longuest part, before the last 50 aa
                                 # OUTPUTs "05_CDS_aa" & "05_CDS_nuc" => NOT INCLUDE THIS CRITERIA
                                 # OUTPUTs "06_CDS_with_M_aa" & "06_CDS_with_M_nuc" => INCLUDE THIS CRITERIA

import string, os, time, re, zipfile, sys, argparse
from dico import dico

def code_universel(F1):
    """ Creates bash for genetic code (key : codon ; value : amino-acid) """
    bash_codeUniversel = {}

    with open(F1, "r") as file:
        for line in file.readlines():
            L1 = string.split(line, " ")
            length1 = len(L1)
            if length1 == 3:
                key = L1[0]
                value = L1[2][:-1]
                bash_codeUniversel[key] = value
            else:
                key = L1[0]
                value = L1[2]
                bash_codeUniversel[key] = value

    return(bash_codeUniversel)

def multiple3(seq):
    """ Tests if the sequence is a multiple of 3, and if not removes extra-bases 
        !! Possible to lost a codon, when I test ORF (as I will decay the ORF) """
    
    m = len(seq)%3
    if m != 0 :
        return seq[:-m], m
    else :
        return seq, m

def detect_Methionine(seq_aa, Ortho, minimal_cds_length):
    """ Detects if methionin in the aa sequence """

    ln = len(seq_aa)
    CUTOFF_Last_50aa = ln - minimal_cds_length

    # Find all indices of occurances of "M" in a string of aa   
    list_indices = [pos for pos, char in enumerate(seq_aa) if char == "M"]

    # If some "M" are present, find whether the first "M" found is not in the 50 last aa (indice < CUTOFF_Last_50aa) ==> in this case: maybenot a CDS
    if list_indices != []:
        first_M = list_indices[0]
        if first_M < CUTOFF_Last_50aa:
            Ortho = 1  # means orthologs found

    return(Ortho)

def ReverseComplement2(seq):
    """ Reverse complement DNA sequence """
    seq1 = 'ATCGN-TAGCN-atcgn-tagcn-'
    seq_dict = { seq1[i]:seq1[i+6] for i in range(24) if i < 6 or 12<=i<=16 }

    return "".join([seq_dict[base] for base in reversed(seq)])

def simply_get_ORF(seq_dna, gen_code):
    seq_by_codons = [seq_dna.upper().replace('T', 'U')[i:i+3] for i in range(0, len(seq_dna), 3)]
    seq_by_aa = [gen_code[codon] if codon in gen_code.keys() else '?' for codon in seq_by_codons]

    return ''.join(seq_by_aa)

def find_good_ORF_criteria_3(bash_aligned_nc_seq, bash_codeUniversel, minimal_cds_length, min_spec):
    # Multiple sequence based : Based on the alignment of several sequences (orthogroup)
    # Criteria 1 : Get the segment in the alignment with no codon stop

    # 1 - Get the list of aligned aa seq for the 3 ORF:
    bash_of_aligned_aa_seq_3ORF = {}
    bash_of_aligned_nuc_seq_3ORF = {}
    BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION = []
    
    for fasta_name in bash_aligned_nc_seq.keys():
        # Get sequence, chek if multiple 3, then get 6 orfs
        sequence_nc = bash_aligned_nc_seq[fasta_name]
        new_sequence_nc, modulo = multiple3(sequence_nc)
        new_sequence_rev = ReverseComplement2(new_sequence_nc)        
        # For each seq of the multialignment => give the 6 ORFs (in nuc)
        bash_of_aligned_nuc_seq_3ORF[fasta_name] = [new_sequence_nc, new_sequence_nc[1:-2], new_sequence_nc[2:-1], new_sequence_rev, new_sequence_rev[1:-2], new_sequence_rev[2:-1]]

        seq_prot_ORF1 = simply_get_ORF(new_sequence_nc, bash_codeUniversel)
        seq_prot_ORF2 = simply_get_ORF(new_sequence_nc[1:-2], bash_codeUniversel)
        seq_prot_ORF3 = simply_get_ORF(new_sequence_nc[2:-1], bash_codeUniversel)
        seq_prot_ORF4 = simply_get_ORF(new_sequence_rev, bash_codeUniversel)
        seq_prot_ORF5 = simply_get_ORF(new_sequence_rev[1:-2], bash_codeUniversel)
        seq_prot_ORF6 = simply_get_ORF(new_sequence_rev[2:-1], bash_codeUniversel)
        
        # For each seq of the multialignment => give the 6 ORFs (in aa)
        bash_of_aligned_aa_seq_3ORF[fasta_name] = [seq_prot_ORF1, seq_prot_ORF2, seq_prot_ORF3, seq_prot_ORF4, seq_prot_ORF5, seq_prot_ORF6]

    # 2 - Test for the best ORF (Get the longuest segment in the alignment with no codon stop ... for each ORF ... the longuest should give the ORF)
    BEST_MAX = 0

    for i in [0,1,2,3,4,5]: # Test the 6 ORFs
        ORF_Aligned_aa = []
        ORF_Aligned_nuc = []

        # 2.1 - Get the alignment of sequence for a given ORF
        # Compare the 1rst ORF between all sequence => list them in ORF_Aligned_aa // them do the same for the second ORF, and them the 3rd
        for fasta_name in bash_of_aligned_aa_seq_3ORF.keys():
            ORFsequence = bash_of_aligned_aa_seq_3ORF[fasta_name][i]
            aa_length = len(ORFsequence)
            ORF_Aligned_aa.append(ORFsequence)   ### List of all sequences in the ORF nb "i" =

        n = i+1

        for fasta_name in bash_of_aligned_nuc_seq_3ORF.keys():
            ORFsequence = bash_of_aligned_nuc_seq_3ORF[fasta_name][i]
            nuc_length = len(ORFsequence)
            ORF_Aligned_nuc.append(ORFsequence) # List of all sequences in the ORF nb "i" =

        # 2.2 - Get the list of sublist of positions whithout codon stop in the alignment
        # For each ORF, now we have the list of sequences available (i.e. THE ALIGNMENT IN A GIVEN ORF)
        # Next step is to get the longuest subsequence whithout stop
        # We will explore the presence of stop "*" in each column of the alignment, and get the positions of the segments between the positions with "*"
        MAX_LENGTH = 0
        LONGUEST_SEGMENT_UNSTOPPED = ""
        j = 0 # Start from first position in alignment
        List_of_List_subsequences = []
        List_positions_subsequence = []
        while j < aa_length:
                column = []
                for seq in ORF_Aligned_aa:
                    column.append(seq[j])
                j = j+1
                if "*" in column:
                    List_of_List_subsequences.append(List_positions_subsequence) # Add previous list of positions
                    List_positions_subsequence = []                              # Re-initialyse list of positions
                else:
                    List_positions_subsequence.append(j)

        # 2.3 - Among all the sublists (separated by column with codon stop "*"), get the longuest one (BETTER SEGMENT for a given ORF)
        LONGUEST_SUBSEQUENCE_LIST_POSITION = []
        MAX=0
        for sublist in List_of_List_subsequences:
            if len(sublist) > MAX and len(sublist) > minimal_cds_length:
                MAX = len(sublist)
                LONGUEST_SUBSEQUENCE_LIST_POSITION = sublist

        # 2.4. - Test if the longuest subsequence start exactly at the beginning of the original sequence (i.e. means the ORF maybe truncated)
        if LONGUEST_SUBSEQUENCE_LIST_POSITION != []:
            if LONGUEST_SUBSEQUENCE_LIST_POSITION[0] == 0:
                CDS_maybe_truncated = 1
            else:
                CDS_maybe_truncated = 0
        else:
            CDS_maybe_truncated = 0


        # 2.5 - Test if this BETTER SEGMENT for a given ORF, is the better than the one for the other ORF (GET THE BEST ORF)
        # Test whether it is the better ORF
        if MAX > BEST_MAX:
            BEST_MAX = MAX
            BEST_ORF = i+1
            BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION = LONGUEST_SUBSEQUENCE_LIST_POSITION


    # 3 - ONCE we have this better segment (BEST CODING SEGMENT)
    # ==> GET THE STARTING and ENDING POSITIONS (in aa position and in nuc position)
    # And get the INDEX of the best ORF [0, 1, or 2]
    if BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION != []:
        pos_MIN_aa = BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION[0]
        pos_MIN_aa = pos_MIN_aa - 1
        pos_MAX_aa = BEST_LONGUEST_SUBSEQUENCE_LIST_POSITION[-1]


        BESTORF_bash_of_aligned_aa_seq = {}
        BESTORF_bash_of_aligned_aa_seq_CODING = {}
        for fasta_name in bash_of_aligned_aa_seq_3ORF.keys():
            index_BEST_ORF = BEST_ORF-1  # cause list going from 0 to 2 in LIST_3_ORF, while the ORF nb is indexed from 1 to 3
            seq = bash_of_aligned_aa_seq_3ORF[fasta_name][index_BEST_ORF]
            seq_coding = seq[pos_MIN_aa:pos_MAX_aa]
            BESTORF_bash_of_aligned_aa_seq[fasta_name] = seq
            BESTORF_bash_of_aligned_aa_seq_CODING[fasta_name] = seq_coding

        # 4 - Get the corresponding position (START/END of BEST CODING SEGMENT) for nucleotides alignment
        pos_MIN_nuc = pos_MIN_aa * 3
        pos_MAX_nuc = pos_MAX_aa * 3

        BESTORF_bash_aligned_nc_seq = {}
        BESTORF_bash_aligned_nc_seq_CODING = {}
        for fasta_name in bash_aligned_nc_seq.keys():
            seq = bash_of_aligned_nuc_seq_3ORF[fasta_name][index_BEST_ORF]
            seq_coding = seq[pos_MIN_nuc:pos_MAX_nuc]
            BESTORF_bash_aligned_nc_seq[fasta_name] = seq
            BESTORF_bash_aligned_nc_seq_CODING[fasta_name] = seq_coding

    else: # no CDS found
        BESTORF_bash_aligned_nc_seq = {}
        BESTORF_bash_aligned_nc_seq_CODING = {}
        BESTORF_bash_of_aligned_aa_seq = {}
        BESTORF_bash_of_aligned_aa_seq_CODING ={}

    # Check whether their is a "M" or not, and if at least 1 "M" is present, that it is not in the last 50 aa
    
    BESTORF_bash_of_aligned_aa_seq_CDS_with_M = {}
    BESTORF_bash_of_aligned_nuc_seq_CDS_with_M = {}

    Ortho = 0
    for fasta_name in BESTORF_bash_of_aligned_aa_seq_CODING.keys():
        seq_aa = BESTORF_bash_of_aligned_aa_seq_CODING[fasta_name]
        Ortho = detect_Methionine(seq_aa, Ortho, minimal_cds_length)   ### DEF6 ###

    # CASE 1: A "M" is present and correctly localized (not in last 50 aa)
    if Ortho == 1:
        BESTORF_bash_of_aligned_aa_seq_CDS_with_M = BESTORF_bash_of_aligned_aa_seq_CODING
        BESTORF_bash_of_aligned_nuc_seq_CDS_with_M = BESTORF_bash_aligned_nc_seq_CODING

    # CASE 2: in case the CDS is truncated, so the "M" is maybe missing:
    if Ortho == 0 and CDS_maybe_truncated == 1:
        BESTORF_bash_of_aligned_aa_seq_CDS_with_M = BESTORF_bash_of_aligned_aa_seq_CODING
        BESTORF_bash_of_aligned_nuc_seq_CDS_with_M = BESTORF_bash_aligned_nc_seq_CODING

    # CASE 3: CDS not truncated AND no "M" found in good position (i.e. before the last 50 aa):
        ## => the 2 bash "CDS_with_M" are left empty ("{}")

    return(BESTORF_bash_aligned_nc_seq,  BESTORF_bash_aligned_nc_seq_CODING, BESTORF_bash_of_aligned_nuc_seq_CDS_with_M, BESTORF_bash_of_aligned_aa_seq, BESTORF_bash_of_aligned_aa_seq_CODING, BESTORF_bash_of_aligned_aa_seq_CDS_with_M)

def write_output_file(results_dict, name_elems, path_out):
    if results_dict != {}:
        name_elems[3] = str(len(results_dict.keys()))
        new_name = "_".join(name_elems)

        out1 = open("%s/%s" %(path_out,new_name), "w")
        for fasta_name in results_dict.keys():
            seq = results_dict[fasta_name]
            out1.write("%s\n" %fasta_name)
            out1.write("%s\n" %seq)
        out1.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("codeUniversel", help="File describing the genetic code (code_universel_modified.txt")
    parser.add_argument("min_cds_len", help="Minmal length of a CDS (in amino-acids)", type=int)    
    parser.add_argument("min_spec", help="Minimal number of species per alignment")
    parser.add_argument("list_files", help="File with all input files names")
    args = parser.parse_args()

    minimal_cds_length = int(args.min_cds_len)  # in aa number
    bash_codeUniversel = code_universel(args.codeUniversel)
    minimum_species = int(args.min_spec)
    
    # Inputs from file containing list of species
    list_files = []
    with open(args.list_files, 'r') as f:
        for line in f.readlines():
            list_files.append(line.strip('\n'))

    # Directories for results
    dirs = ["04_BEST_ORF_nuc", "04_BEST_ORF_aa", "05_CDS_nuc", "05_CDS_aa", "06_CDS_with_M_nuc", "06_CDS_with_M_aa"]
    for directory in dirs:
        os.mkdir(directory)

    count_file_processed, count_file_with_CDS, count_file_without_CDS, count_file_with_CDS_plus_M = 0, 0, 0, 0
    count_file_with_cds_and_enought_species, count_file_with_cds_M_and_enought_species = 0, 0

    # ! : Currently, files are named "Orthogroup_x_y_sequences.fasta, where x is the number of the orthogroup (not important, juste here to make a distinct name),
    # and y is the number of sequences/species in the group. These files are outputs of blastalign, where species can be removed. y is then modified.
    name_elems = ["orthogroup", "0", "with", "0", "species.fasta"]

    # by fixing the counter here, there will be some "holes" in the outputs directories (missing numbers), but the groups between directories will correspond
    #n0 = 0

    for file in list_files:
        #n0 += 1

        count_file_processed = count_file_processed + 1
        nb_gp = file.split('_')[1] # Keep trace of the orthogroup number
        fasta_file_path = "./%s" %file    
        bash_fasta = dico(fasta_file_path)
        BESTORF_nuc, BESTORF_nuc_CODING, BESTORF_nuc_CDS_with_M, BESTORF_aa, BESTORF_aa_CODING, BESTORF_aa_CDS_with_M  = find_good_ORF_criteria_3(bash_fasta, bash_codeUniversel, minimal_cds_length, minimum_species)
        
        name_elems[1] = nb_gp

        # Update counts and write group in corresponding output directory
        if BESTORF_nuc != {}: 
            count_file_with_CDS += 1
            if len(BESTORF_nuc.keys()) >= minimum_species :
                count_file_with_cds_and_enought_species += 1
                write_output_file(BESTORF_nuc, name_elems, dirs[0]) # OUTPUT BESTORF_nuc
                write_output_file(BESTORF_aa, name_elems, dirs[1]) # The most interesting
        else:
            count_file_without_CDS += 1

        if BESTORF_nuc_CODING != {} and len(BESTORF_nuc_CODING.keys()) >= minimum_species:
            write_output_file(BESTORF_nuc_CODING, name_elems, dirs[2])
            write_output_file(BESTORF_aa_CODING, name_elems, dirs[3])

        if BESTORF_nuc_CDS_with_M != {}:
            count_file_with_CDS_plus_M += 1
            if len(BESTORF_nuc_CDS_with_M.keys()) >= minimum_species :
                count_file_with_cds_M_and_enought_species += 1
                write_output_file(BESTORF_nuc_CDS_with_M, name_elems, dirs[4])
                write_output_file(BESTORF_aa_CDS_with_M, name_elems, dirs[5])

    print "*************** CDS detection ***************"
    print "\nFiles processed: %d" %count_file_processed
    print "\tFiles with CDS: %d" %count_file_with_CDS
    print "\tFiles wth CDS and more than %s species: %d" %(minimum_species, count_file_with_cds_and_enought_species)
    print "\t\tFiles with CDS plus M (codon start): %d" %count_file_with_CDS_plus_M
    print "\t\tFiles with CDS plus M (codon start) and more than %s species: %d" %(minimum_species,count_file_with_cds_M_and_enought_species) 
    print "\tFiles without CDS: %d \n" %count_file_without_CDS
    print ""

if __name__ == '__main__':
    main()
