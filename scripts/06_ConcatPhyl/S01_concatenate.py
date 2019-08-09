#!/usr/bin/python
## Author: Eric Fontanillas
## Last modification: 07/2019
## Subject: find and remove indels

import string, os, time, re, sys, itertools

"""
Parse fasta alignment file (with indels) and replace genes headers with matching
species name

Args:
    - file_in (str): fasta file
    - species_id_list (list) : all the species identifier (basically the origin
    name of the fasta file, without .extension)

Return:
    - alignmetn as a dictionary (species: sequence)
"""
def dico(file_in, species_id_list):
    dic = {}

    with open(file_in, "r") as f_in:

        for name, query in itertools.izip_longest(*[f_in]*2):
            # Match gene header with species ID (maybe better with itertools ?)
            for s in species_id_list:
                if re.match('.+'+s, name.rstrip()):
                    dic[s] = query.rstrip()

    return dic

"""
Parses all alignment files and concatenate all sequences to make a 'super-alignment'

Args:
    - list_files (list) all the files names (retrived previously from an input file)
    - species_id_list (list) : all the species identifier (basically the origin
    name of the fasta file, without .extension)
"""
def concatenate(list_files, species_id_list):

    concatenated_alignments = {}

    for species_ID in species_id_list:
        concatenated_alignments[species_ID] = ''

    len_concat = 0
    nb_locus = 0
    pos = 1
    list_genes_position = []

    for file in list_files:
        nb_locus += 1

        # Open alignments and make headers match with files names
        dico_seq = dico(file, species_id_list)

        # Get alignment length (file and total)
        ln = len(dico_seq.values()[0])
        len_concat = len_concat + ln

        # Get genes positions for RAxML
        pos_start = pos
        pos_end = pos + ln - 1
        pos = pos_end + 1 # update pos for next file
        position = "%d-%d" %(pos_start, pos_end)
        sublist = [file, position]
        list_genes_position.append(sublist)

        # Get missing species in the alignment
        for sp in species_id_list:
            if sp not in dico_seq.keys():
                dico_seq[sp] = '-' * ln

        # Concatenate
        for k,v in dico_seq.items():
            concatenated_alignments[k] += v

    return(concatenated_alignments, len_concat, nb_locus, list_genes_position)


""" get codon position """
def get_codon_position(seq_inORF):

    ln = len(seq_inORF)

    i=0
    seq_pos1=""
    seq_pos2=""
    seq_pos12=""
    seq_pos3=""
    while i<ln:
       pos1 =  seq_inORF[i]
       pos2 =  seq_inORF[i+1]
       pos3 =  seq_inORF[i+2]

       seq_pos1 = seq_pos1 + pos1
       seq_pos2 = seq_pos2 + pos2
       seq_pos12 = seq_pos12 + pos1 + pos2
       seq_pos3 = seq_pos3 + pos3

       i = i+3

    return(seq_pos1, seq_pos2, seq_pos12, seq_pos3)

#######################
##### RUN RUN RUN #####
#######################

list_species = []
species_id_list = []
fasta = "^.*fasta$"
i=3

## Arguments
infiles_filter_assemblies = sys.argv[1]
format_run = sys.argv[2]

## add file to list_species
list_species = str.split(infiles_filter_assemblies,",")

# In species_id_list, record the file name trimmed of the extension,
# to match gene identifiers, ids common part being the oriinal file name
# -> FileNames MUST NOT have a '.' other than the one defining the file extension
for name in list_species :
    name = name.split('.')[0]
    species_id_list.append(name)

# add alignment files
list_files = []
with open(sys.argv[3], 'r') as f:
    for line in f.readlines():
        list_files.append(line.rstrip())


### 1 ### Proteic
if format_run == "proteic" :

    OUT1 = open("concatenation_aa.fasta", "w")
    OUT2 = open("concatenation_aa.phy", "w")
    OUT3 = open("concatenation_aa.nxs", "w")
    OUT_PARTITION_gene_AA = open("partitions_gene_AA","w")

    # Get bash with concatenation
    bash_concatenation, ln, nb_locus, list_genes_position= concatenate(list_files, species_id_list)

    ##Write gene AA partition file for RAxML
    for sublist in list_genes_position:
        name = sublist[0]
        positions=sublist[1]
        OUT_PARTITION_gene_AA.write("DNA,%s=%s\n"%(name,positions))
    OUT_PARTITION_gene_AA.close()

    # Get "ntax" for NEXUS HEADER
    nb_taxa = len(bash_concatenation.keys())

    print "******************** CONCATENATION ********************\n"
    print "Process amino-acid concatenation:"
    print "\tNumber of taxa aligned = %d" %nb_taxa
    print "\tNumber of loci concatenated = %d\n" %nb_locus
    print "\tTotal length of the concatenated sequences = %d" %ln

    # nexus header
    OUT3.write("#NEXUS\n\n")
    OUT3.write("Begin data;\n")
    OUT3.write("\tDimensions ntax=%d nchar=%d;\n" %(nb_taxa, ln))
    OUT3.write("\tFormat datatype=aa gap=-;\n")
    OUT3.write("\tMatrix\n")

    # phylip header
    OUT2.write("   %d %d\n" %(nb_taxa, ln))

    # write outputs
    for seq_name in bash_concatenation.keys():
        seq = bash_concatenation[seq_name]

        # Filtering the sequence in case of remaining "?"
        seq = string.replace(seq, "?", "-")

        # fasta format
        OUT1.write(">%s\n" %seq_name)
        OUT1.write("%s\n" %seq)

        # phylip format
        OUT2.write("%s\n" %seq_name)
        OUT2.write("%s\n" %seq)

        # nexus
        OUT3.write("%s" %seq_name)
        OUT3.write("      %s\n" %seq)

    OUT3.write("\t;\n")
    OUT3.write("End;\n")
    OUT1.close()
    OUT2.close()
    OUT2.close()

### 2 ### Nucleic
elif format_run == "nucleic" :

    OUT1 = open("concatenation_nuc.fasta", "w")
    OUT2 = open("concatenation_nuc.phy", "w")
    OUT3 = open("concatenation_nuc.nxs", "w")

    OUT1_pos12 = open("concatenation_pos12_nuc.fasta", "w")
    OUT2_pos12 = open("concatenation_pos12_nuc.phy", "w")
    OUT3_pos12 = open("concatenation_pos12_nuc.nxs", "w")

    OUT1_pos3 = open("concatenation_pos3_nuc.fasta", "w")
    OUT2_pos3 = open("concatenation_pos3_nuc.phy", "w")
    OUT3_pos3 = open("concatenation_pos3_nuc.nxs", "w")

    OUT_PARTITION_codon_12_3 = open("partitions_codon12_3","w")
    OUT_PARTITION_gene_NUC = open("partitions_gene_NUC","w")
    OUT_PARTITION_gene_PLUS_codon_12_3 = open("partitions_gene_PLUS_codon12_3","w")

    ## Get bash with concatenation
    bash_concatenation, ln, nb_locus, list_genes_position = concatenate(list_files, species_id_list)

    # for k,v in bash_concatenation.items():
        # print k,v

    ln_12 = ln/3*2   ### length of the alignment when only the 2 first codon position
    ln_3 = ln/3      ### length of the alignment when only the third codon position

    ## Write partition files for RAxML subsequent runs
    # a # Codon partition
    OUT_PARTITION_codon_12_3.write("DNA, p1=1-%d\\3,2-%d\\3\n" %(ln, ln))
    OUT_PARTITION_codon_12_3.write("DNA, p2=3-%d\\3\n" %(ln))
    OUT_PARTITION_codon_12_3.close()

    # b # Gene partition
    for sublist in list_genes_position:
        name=sublist[0]
        positions=sublist[1]
        OUT_PARTITION_gene_NUC.write("DNA,%s=%s\n"%(name,positions))
    OUT_PARTITION_gene_NUC.close()

    # c # Mixed partition (codon + gene)
    for sublist in list_genes_position:
        name = sublist[0]
        positions = sublist[1]
        S1 = string.split(positions, "-")
        pos_start1 = string.atoi(S1[0])
        pos_end = string.atoi(S1[1])
        pos_start2=pos_start1+1
        pos_start3=pos_start2+1
        partition1 = "DNA, %s_1=%d-%d\\3,%d-%d\\3\n" %(name,pos_start1, pos_end, pos_start2, pos_end)
        partition2 = "DNA, %s_2=%d-%d\\3\n" %(name,pos_start3, pos_end)
        OUT_PARTITION_gene_PLUS_codon_12_3.write(partition1)
        OUT_PARTITION_gene_PLUS_codon_12_3.write(partition2)

    OUT_PARTITION_gene_PLUS_codon_12_3.close()

    ## Get "ntax" for NEXUS HEADER
    nb_taxa = len(bash_concatenation.keys())

    print "******************** CONCATENATION ********************\n"
    print "Process nucleotides concatenation:"
    print "\tNumber of taxa aligned = %d" %nb_taxa
    print "\tNumber of loci concatenated = %d\n" %nb_locus
    print "\tTotal length of the concatenated sequences [All codon positions] = %d" %ln
    print "\t\tTotal length of the concatenated sequences [Codon positions 1 & 2] = %d" %ln_12
    print "\t\tTotal length of the concatenated sequences [Codon position 3] = %d" %ln_3

    # nexus header
    OUT3.write("#NEXUS\n\n")
    OUT3.write("Begin data;\n")
    OUT3.write("\tDimensions ntax=%d nchar=%d;\n" %(nb_taxa, ln))
    OUT3.write("\tFormat datatype=dna gap=-;\n")
    OUT3.write("\tMatrix\n")

    OUT3_pos12.write("#NEXUS\n\n")
    OUT3_pos12.write("Begin data;\n")
    OUT3_pos12.write("\tDimensions ntax=%d nchar=%d;\n" %(nb_taxa, ln_12))
    OUT3_pos12.write("\tFormat datatype=dna gap=-;\n")
    OUT3_pos12.write("\tMatrix\n")

    OUT3_pos3.write("#NEXUS\n\n")
    OUT3_pos3.write("Begin data;\n")
    OUT3_pos3.write("\tDimensions ntax=%d nchar=%d;\n" %(nb_taxa, ln_3))
    OUT3_pos3.write("\tFormat datatype=dna gap=-;\n")
    OUT3_pos3.write("\tMatrix\n")

    # phylip header
    OUT2.write("   %d %d\n" %(nb_taxa, ln))
    OUT2_pos12.write("   %d %d\n" %(nb_taxa, ln_12))
    OUT2_pos3.write("   %d %d\n" %(nb_taxa, ln_3))

    ## Print outputs
    for seq_name in bash_concatenation.keys():
        seq = bash_concatenation[seq_name]

        # Filtering the sequence in case of remaining "?"
        seq = string.replace(seq, "?", "-")

        # Get the differentes codons partitions
        seq_pos1, seq_pos2, seq_pos12, seq_pos3 = get_codon_position(seq)

        # fasta
        OUT1.write(">%s\n" %seq_name)
        OUT1.write("%s\n" %seq)
        OUT1_pos12.write(">%s\n" %seq_name)
        OUT1_pos12.write("%s\n" %seq_pos12)
        OUT1_pos3.write(">%s\n" %seq_name)
        OUT1_pos3.write("%s\n" %seq_pos3)

        # phylip
        OUT2.write("%s\n" %seq_name)
        OUT2.write("%s\n" %seq)
        OUT2_pos12.write("%s\n" %seq_name)
        OUT2_pos12.write("%s\n" %seq_pos12)
        OUT2_pos3.write("%s\n" %seq_name)
        OUT2_pos3.write("%s\n" %seq_pos3)

        # nexus
        OUT3.write("%s" %seq_name)
        OUT3.write("      %s\n" %seq)
        OUT3_pos12.write("%s" %seq_name)
        OUT3_pos12.write("      %s\n" %seq_pos12)
        OUT3_pos3.write("%s" %seq_name)
        OUT3_pos3.write("      %s\n" %seq_pos3)

    OUT3.write("\t;\n")
    OUT3.write("End;\n")
    OUT3_pos12.write("\t;\n")
    OUT3_pos12.write("End;\n")
    OUT3_pos3.write("\t;\n")
    OUT3_pos3.write("End;\n")

    OUT1.close()
    OUT2.close()
    OUT3.close()
    OUT1_pos12.close()
    OUT2_pos12.close()
    OUT3_pos12.close()
    OUT1_pos3.close()
    OUT2_pos3.close()
    OUT3_pos3.close()

print "\n\n\n******************** RAxML RUN ********************\n"
