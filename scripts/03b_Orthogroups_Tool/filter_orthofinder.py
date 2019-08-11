#!/usr/bin/env python
# Commandline : ./filter_orthofinder.py <Orthogroups.txt> <Nb_of_studied_species> <minimal_nb_species_per_group> [-v] [-p]

## This script takes an output file of OrthoFinder (Orthogroups.txt), which contains a set of orthogroups,
## and rewrite it to split each orthogroup into a single fasta file.

import os, string, glob, argparse, csv, itertools, re
import numpy as np
import pandas as pd

## PRELIMINARY :
""" rewrite TransDecoder headers within Orthogroups.txt """
def rewriteHeaders(path):
    out = open('tmp_orthogroups.txt', 'w')

    with open(path, 'r') as f:
        for line in f.readlines():
            l = line.rstrip()
            group = re.findall(r"Gene\.[0-9]+__(.*?)__g\.[0-9]+", l)#[::2]

            for g in group:
                out.write(g + ' ')
            out.write('\n')

    out.close()

## PART 1 : Make a dictionary of {IDs : sequence}

""" Build a hash table with gene IDs and gene sequences from fasta made from input files
    Returns a dictionnary """
def hashSequences(path):
    hashTable = {}
    # WARNING : sequences are expected to be on one line. If not, biopython can do it
    for file in path:
        with open(file, "r") as origin:
            for line1,line2 in itertools.izip_longest(*[origin]*2):
                gene = line1.rstrip().strip('>')
                sequence = line2.rstrip().strip('>')
                hashTable[gene] = sequence

    return hashTable

## PART 2 : Create orthogroups file (one file per orthogroup)

""" Takes a file.txt of orthogroups as parameter and keeps the orthogroups where there are at least args.minspec loci
    WARNING : sequences names within the groups must the same as IDs in fasta files from Filter_Assemblies. if not, the
    dictionnary will be false. That's is why the script "format_transdecoder_headers is for.
    Returns a 2D list (each list being a list of loci) """
def formatAndFilter(orthogroups, mini, nbspecs, hashTable, verbose, paralogs, list_specs):

    """ Builds a 2D array for a summary
        Returns a numpy 2D array """
    def countings(listOrthogroups, nbcols, list_specs):
        #listOrthogroups.sort().reverse()
        #nblines = len(listOrthogroups[0])
        nblines = 0
        for group in listOrthogroups:
            if len(group) > nblines:
                nblines = len(group)
        matrix = np.array([[0]*nbcols]*nblines)
        # empty lines are avoided : first line of the frame is the line for minimal number of sequences in a group (>=mini)
        # for now, this feature diseappear when using numpy arrays and pandas :/

        for group in listOrthogroups:
            listSpecs = []
            for loci in group:
                for spec in list_specs:
                    if re.match('>'+spec, loci) and spec not in listSpecs:
                        listSpecs.append(spec)
                        break
                # if loci.getHeader()[1:3] not in listSpecs:
                #     listSpecs.append(loci.getHeader()[1:3])
            matrix[len(group)-1][len(listSpecs)-1] += 1

        return matrix

    """ numpy 2D array in a nice dataframe
        Returns a pandas 2D dataframe """
    def asFrame(matrix) :
        index = [0]*len(matrix)
        colnames = [0]*len(matrix[0])
        index = [str(i+1)+" seqs" for i in range(len(matrix))]
        colnames = [str(i+1)+" sps" for i in range(len(matrix[0]))]
        df = pd.DataFrame(matrix, index=index, columns=colnames)
        return df # Mettre une selection pour ne renvoyer que les lignes et les colonnes qui somment > 0
        #return df.loc['4 seqs':'9 seqs'].loc[:,colnames[3:]]

    """ Writes each orthogroup in a fasta file. Retrieves sequences with a hash table """
    def writeOutputFile(orthogroup, hashTable, i, naming):
        length = len(orthogroup)
        name =""
        if naming:
            name="orthogroup_{}_{}_sequences_withParalogs.fasta".format(i, length)
        else :
            name = "orthogroup_{}_{}_sequences.fasta".format(i, length)
        result = open(name, "w")
        with result:
            for locus in orthogroup:
                result.write(">{}\n".format(locus)) # write geneID. ">%s\n" before
                result.write("{}\n".format(hashTable[locus])) # write sequence

    # FUNCTION

    orthogroups = open(orthogroups, "r")
    list_orthogroups = []
    # STEP 1 - Read file into a list ----------------------------------------------
    with orthogroups:
        while (1):
            group = orthogroups.readline().rstrip()
            if not group:
                break
            else:
                list_orthogroups.append(group)

    # STEP 2 - Paralogous filtering -----------------------------------------------
    if verbose or paralogs:
        list_orthogroups_withpara = []
    list_orthogroups_format = []

    i,j = 1,1
    for group in list_orthogroups:
        group = string.split(group, " ") # list of lists
        group.sort()
        if verbose and len(group) >= mini:
            list_orthogroups_withpara.append(group)
        if paralogs and len(group) >= mini:
            writeOutputFile(group, hashTable, j, True)
        j += 1
        new_group = []
        rang=-1
        # Keep only one paralogs per species (1st encounter)
        species = []
        for loci in group:
            # detect species
            sp = ''
            for spec in list_specs:
                if re.match(spec, loci):
                    sp = spec
                    break
            # record loci if species not already detected
            if sp not in species:
                new_group.append(loci)
                species.append(sp)

        if len(new_group) >= mini: # Drop too small orthogroups
            list_orthogroups_format.append(new_group)
            writeOutputFile(new_group, hashTable, i, False)
            i += 1

    # STEP 3 - Print summaries ----------------------------------------------------
    if verbose:
        print "  Summary before paralogous filtering : \n"
        df1 = asFrame(countings(list_orthogroups_withpara, nbspecs, list_specs))
        print df1.loc[df1.ne(0).any(1),df1.ne(0).any()]
        #print "  Summary before paralogous filtering : \n",countings(list_orthogroups_withpara, nbspecs),"\n"
    print "\n  Summary after paralogous filtering : \n"
    df2 = asFrame(countings(list_orthogroups_format, nbspecs, list_specs))
    print df2.loc[df2.ne(0).any(1),df2.ne(0).any()]

    return len(list_orthogroups_format) #list_orthogroups_no_para

## MAIN

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("orthogroups", help="Orthogroups.txt file from OrthoFinder")
    parser.add_argument("infiles", help="fasta files used for OrthoFinder (after or before TransDecoder), separated by commas")
    parser.add_argument("list_specs", help="A text file with the list of all species names (i.e original files names from FIlter Assemblies, without the '.fasta'")
    parser.add_argument("minspec", type=int, help="Minimal number of species to keep per group")
    parser.add_argument("-v", "--verbose", action="store_true", help="Add another summary table : countings before paralogous genes filtering")
    parser.add_argument("-p", "--paralogs", action="store_true", help="Proceeds to write orthogroups also before paralogous filtering")
    args = parser.parse_args()

    print "\n-This script works on the 'Orthogroups' file output of Orthofinder to split each orthogroup in a single fasta file."
    print "-It also gets rid of orthogroups with less sequences than the number specified by the user."

    # Headers
    print "Re-writing Orthogroups headers ...\n"
    rewriteHeaders(args.orthogroups)

    # Build hashtable
    print "  Building hashTable IDs/sequences ...\n"
    infiles_good = str.split(args.infiles, ",")
    nbspec = len(infiles_good)
    #path = glob.glob('*.fasta')
    hashTable = hashSequences(infiles_good)

    list_specs = []
    with open(args.list_specs) as ls:
        for line in ls.readlines():
            list_specs.append(line.rstrip())

    # Open txt file with orthogroups
    print "  Reading Orthogroups.txt and writing orthogroups to separated files..."
    print "    (Dropping orthogroups of less than {} loci.)\n".format(args.minspec)
    list_orthogroups = formatAndFilter('tmp_orthogroups.txt', args.minspec, nbspec, hashTable, args.verbose, args.paralogs, list_specs)
    print "\n{} filtered orthogroups have been written in separated files".format(list_orthogroups)

    # Move output files in a new directory
    if args.paralogs:
        os.system("mkdir orthogroups_withParalogs")
        path=glob.glob("*_withParalogs.fasta")
        for file in path:
            os.system("mv {} orthogroups_withParalogs".format(file))

    os.system("mkdir filtered_orthogroups")
    path = glob.glob("*_sequences.fasta")
    for file in path:
        os.system("mv {} filtered_orthogroups".format(file))

    print "  \nFiltered orthogroups are written in the directory 'filtered_orthogroups'"
    if args.paralogs:
        print "  \nFull orthogroups files are written in the directory 'orthogroups_withParalogs'\n"

if __name__ == "__main__":
    main()
