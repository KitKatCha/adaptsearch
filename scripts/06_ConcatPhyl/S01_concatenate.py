#!/usr/bin/python
## Author: Eric Fontanillas
## Last modification: 17/06/2011
## Subject: find and remove indels

###############################
##### DEF 0 : Dico fasta  #####
###############################
def dico(F2):
    dicoco = {}
    with open(F2, "r") as file:
        for name, query in itertools.izip_longest(*[file]*2):
            if not name:
                break
            if name[0] == ">":
                fasta_name_query = name[:-1]
                Sn = string.split(fasta_name_query, "||")
                fasta_name_query = Sn[0]
                fasta_seq_query = query[:-1]
                dicoco[fasta_name_query]=fasta_seq_query
    return dicoco
###################################################################################


####################
###### DEF 11 ######
####################
## Concatenate sequences
###########################
def concatenate(L_IN, SPECIES_ID_LIST):
    ## 4 ## Process files
    ## 4.1 ## Create the bash and the fasta names entries (name of the species)
    bash_concat = {}

    for species_ID in SPECIES_ID_LIST:
        bash_concat[species_ID] = ''

    ln_concat = 0
    nb_locus = 0
    pos=1
    list_genes_position=[]
    ## 4.2 ## Concatenate
    for file in L_IN:            
        if os.stat(file).st_size == 0:            
            pass
        else :
            nb_locus=nb_locus+1

            ## a ## Open alignments        
            dico_seq = dico(file)   ### DEF 0 ###        
            ## b ## Get alignment length + genes positions for RAxML
            key0 = dico_seq.keys()[0]
            ln = len(dico_seq[key0])
            ln_concat = ln_concat + ln

            pos_start = pos
            pos_end = pos+ln-1
            pos=pos_end+1
            position="%d-%d" %(pos_start, pos_end)
            RAxML_name = file[:-6]
            sublist = [RAxML_name, position]
            list_genes_position.append(sublist)

            ## c ## Generate "empty" sequence with alignment length * "-"
            empty_seq = "-" * ln

            ## d ## Concatenate
            ## d.1 ## Detect missing species in this alignment
            list_ID=[]
            list_absent_ID=[]
            bash_fastaName={}
            for fasta_name in dico_seq:
                ID = fasta_name[1:3]
                list_ID.append(ID)
                seq = dico_seq[fasta_name]
                bash_fastaName[ID]=fasta_name
            for sp_ID in SPECIES_ID_LIST:
                if sp_ID not in list_ID:
                    list_absent_ID.append(sp_ID)

            for ID in SPECIES_ID_LIST:
                if ID in list_absent_ID:
                    bash_concat[ID] = bash_concat[ID] + empty_seq
                else:
                    fasta_name = bash_fastaName[ID]
                    seq = dico_seq[fasta_name]
                    bash_concat[ID] = bash_concat[ID] + seq

    return(bash_concat, ln_concat, nb_locus, list_genes_position)
####################################


########################################
##### DEF 12 : get codon position  #####
########################################
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
###############################################################################



#######################
##### RUN RUN RUN #####
#######################
import string, os, time, re, sys, itertools

list_species = []
SPECIES_ID_LIST = []
fasta = "^.*fasta$"
i=3

## Arguments
infiles_filter_assemblies = sys.argv[1]
format_run = sys.argv[2]
input_alignments = sys.argv[3]

## add file to list_species
list_species = str.split(infiles_filter_assemblies,",")

## in SPECIES_ID_LIST, only the 2 first letters of name of species
for name in list_species :
    name = name[:2]
    SPECIES_ID_LIST.append(name)

## add alignment files to L_IN
L_IN = str.split(input_alignments,",")
print(L_IN)

### 1 ### Proteic
if format_run == "proteic" :

    OUT1 = open("02_Concatenation_aa.fas", "w")
    OUT2 = open("02_Concatenation_aa.phy", "w")
    OUT3 = open("02_Concatenation_aa.nex", "w")
    OUT_PARTITION_gene_AA = open("06_partitions_gene_AA","w")


    ##  Get bash with concatenation
    bash_concatenation, ln, nb_locus,list_genes_position= concatenate(L_IN, SPECIES_ID_LIST)    ### DEF 11 ##

    ## Write gene AA partition file for RAxML
    for sublist in list_genes_position:
        name = sublist[0]
        positions=sublist[1]
        OUT_PARTITION_gene_AA.write("DNA,%s=%s\n"%(name,positions))
    OUT_PARTITION_gene_AA.close()

    ## Get "ntax" for NEXUS HEADER
    nb_taxa = len(bash_concatenation.keys())

    print "******************** CONCATENATION ********************\n"
    print "Process amino-acid concatenation:"
    print "\tNumber of taxa aligned = %d" %nb_taxa
    print "\tNumber of loci concatenated = %d\n" %nb_locus
    print "\tTotal length of the concatenated sequences = %d" %ln


    ## Print NEXUS HEADER:
    OUT3.write("#NEXUS\n\n")
    OUT3.write("Begin data;\n")
    OUT3.write("\tDimensions ntax=%d nchar=%d;\n" %(nb_taxa, ln))
    OUT3.write("\tFormat datatype=aa gap=-;\n")
    OUT3.write("\tMatrix\n")

    ## Print PHYLIP HEADER:
    OUT2.write("   %d %d\n" %(nb_taxa, ln))

    ## 3.5 ## Print outputs
    for seq_name in bash_concatenation.keys():
        seq = bash_concatenation[seq_name]

        ## Filtering the sequence in case of remaining "?"
        seq = string.replace(seq, "?", "-")

        #print seq FASTA FORMAT
        OUT1.write(">%s\n" %seq_name)
        OUT1.write("%s\n" %seq)

        #print seq PHYLIP FORMAT
        OUT2.write("%s\n" %seq_name)
        OUT2.write("%s\n" %seq)

        #print seq NEXUS FORMAT
        OUT3.write("%s" %seq_name)
        OUT3.write("      %s\n" %seq)

    OUT3.write("\t;\n")
    OUT3.write("End;\n")
    OUT1.close()
    OUT2.close()
    OUT2.close()


### 2 ### Nucleic
elif format_run == "nucleic" :

    OUT1 = open("03_Concatenation_nuc.fas", "w")
    OUT2 = open("03_Concatenation_nuc.phy", "w")
    OUT3 = open("03_Concatenation_nuc.nex", "w")

    OUT1_pos12 = open("03_Concatenation_pos12_nuc.fas", "w")
    OUT2_pos12 = open("03_Concatenation_pos12_nuc.phy", "w")
    OUT3_pos12 = open("03_Concatenation_pos12_nuc.nex", "w")

    OUT1_pos3 = open("03_Concatenation_pos3_nuc.fas", "w")
    OUT2_pos3 = open("03_Concatenation_pos3_nuc.phy", "w")
    OUT3_pos3 = open("03_Concatenation_pos3_nuc.nex", "w")

    OUT_PARTITION_codon_12_3 = open("05_partitions_codon12_3","w")
    OUT_PARTITION_gene_NUC = open("05_partitions_gene_NUC","w")
    OUT_PARTITION_gene_PLUS_codon_12_3 = open("05_partitions_gene_PLUS_codon12_3","w")

    ## Get bash with concatenation
    bash_concatenation, ln, nb_locus, list_genes_position = concatenate(L_IN, SPECIES_ID_LIST)    ### DEF 11 ##
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


    ## Print NEXUS HEADER:
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

    ## Print PHYLIP HEADER:
    OUT2.write("   %d %d\n" %(nb_taxa, ln))
    OUT2_pos12.write("   %d %d\n" %(nb_taxa, ln_12))
    OUT2_pos3.write("   %d %d\n" %(nb_taxa, ln_3))

    ## Print outputs
    for seq_name in bash_concatenation.keys():
        seq = bash_concatenation[seq_name]

        ## Filtering the sequence in case of remaining "?"
        seq = string.replace(seq, "?", "-")

        ## Get the differentes codons partitions
        seq_pos1, seq_pos2, seq_pos12, seq_pos3 = get_codon_position(seq)    ### DEF 12 ###

        #print seq FASTA FORMAT
        OUT1.write(">%s\n" %seq_name)
        OUT1.write("%s\n" %seq)
        OUT1_pos12.write(">%s\n" %seq_name)
        OUT1_pos12.write("%s\n" %seq_pos12)
        OUT1_pos3.write(">%s\n" %seq_name)
        OUT1_pos3.write("%s\n" %seq_pos3)

        #print seq PHYLIP FORMAT
        OUT2.write("%s\n" %seq_name)
        OUT2.write("%s\n" %seq)
        OUT2_pos12.write("%s\n" %seq_name)
        OUT2_pos12.write("%s\n" %seq_pos12)
        OUT2_pos3.write("%s\n" %seq_name)
        OUT2_pos3.write("%s\n" %seq_pos3)

        #print seq NEXUS FORMAT
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
