#!/usr/bin/env python
# coding: utf8
# September 2017 - Author : Victor Mataigne (Station Biologique de Roscoff - ABiMS)
# Command line : ./pogsPOO.py <list_of_input_files_separated_by_commas> <minimal number of species per group> [-v) [-p]

"""
What it does:
    - pogs.py parses output files from the "pairwise" tool of the AdaptSearch suite and proceeds to gather genes in orthogroups (using transitivity).
    - A minimal number of species per group has to be set.

Future improvements:
    - Avoid multiple parsing of the orthogroup list and reduce computation time - it is parsed several times :
        - A first time in filterParalogs
        - A second ant a third time in countings
        - Reduce also memory consumption :
            -> the script build several orthogroups lists; if it was replaced by integers or list of integers, or dictionnaries containing countings,
            it would be a better memory usage
"""

import os, argparse
import numpy as np
import pandas as pd

""" Definition of a locus : header + sequence + a tag """
class Locus:    

    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence
        self.tagged = False

    def __str__(self):
        return "{}{}".format(self.header, self.sequence)

    def __eq__(self, other):
        return self.getHeader() == other.getHeader()  # Test if two loci are the same

    def __hash__(self):
        return hash((self.header, self.sequence))  # Make the object iterable and hashable

    def getHeader(self):
        return self.header

    def getSequence(self):
        return self.sequence

    def getTag(self):
        return self.tagged

    def prettyPrint(self): 
        # Used for debugging : print "{ Header : ", self.header[0:-1], "Tag : ", self.tagged, " }"
        print "[ Header : {header} ]".format(header=self.header[0:-1])

    def prettyPrint2(self):
        print "[ Header : {header} Sequence : {sequence} ]".format(header=self.header[0:-1], sequence=self.sequence[0:-1])

""" Applies the getPairwiseCouple() function to a list of files and return a big list with ALL pairwises couples 
    Returns a list of sets (2 items per set) """
def getListPairwiseAll(listPairwiseFiles):    
    
    # Sub-Function

    """ Reads an output file from the 'Pairwise' tool (AdaptSearch suite) and returns its content into a list 
        Returns a list of sets (2 items per set) """
    def getPairwiseCouple(pairwiseFile):        
        list_pairwises_2sp = []
        with open(pairwiseFile, "r") as file:
            while (1):  # Ugly !
                name, sequence, name2, sequence2 = file.readline(), file.readline(), file.readline(), file.readline()
                if not name: break # Use assert ?
                # One locus every two lines (one pairwise couple = 4 lines) : header + sequence
                locus1 = Locus(name, sequence)
                locus2 = Locus(name2, sequence2)
                group = set([])
                group.add(locus1)
                group.add(locus2)
                list_pairwises_2sp.append(group)
        return (list_pairwises_2sp)

    # Function
    list_pairwises_allsp = []
    for file in listPairwiseFiles:
        listPairwises = getPairwiseCouple(file)
        for pairwise in listPairwises:
            list_pairwises_allsp.append(pairwise)  # all pairwises in the same 1D list
    return list_pairwises_allsp

""" Proceeds to create orthogroups by putting together pairwise couples sharing a locus.
    Iterates over the orthogroups list and tag to 'True' the pairwise couple already gathered in a group to avoid
    redondancy. Writes each orthogroup in a fasta file
    Returns an integer (a list length) """
def makeOrthogroups(list_pairwises_allsp, minspec, nb_rbh, verbose):

    # Sub-funtions

    """ Check if a locus/group has already been treated in makeOrthogroups()
        Returns a boolean """
    def checkIfTagged(pair):
        tag = True
        for element in pair:
            if not element.getTag() and tag: # use a list comprehension maybe ?
                tag = False
        return tag

    """ True means a locus/group has already been treated in makeOrthogroups()
        A stronger code would be to implement a method inside the class Locus """
    def tagGroup(pair):
        for element in pair:
            element.tagged = True

    """ Write an orthogroup in a file 
        The previous version had a supplmental boolean parameter, which was used to know if the
        written group was already or not filtered of its paralogs (the 'name' was then modified, 
        as well as the output)"""
    def writeOutputFile(orthogroup, number):
        name = "orthogroup_{}_with_{}_sequences.fasta".format(number, len(orthogroup))
        result = open(name, "w")
        with result:
            for locus in orthogroup:
                if locus.getHeader()[-1] == "\n":
                    result.write("%s" % locus.getHeader())  # write geneID
                else :
                    result.write("%s\n" % locus.Header())  # write geneID
                if locus.getSequence()[-1] == "\n":
                    result.write("%s" % locus.getSequence())  # write sequence
                else :
                    result.write("%s\n" % locus.getSequence())  # write sequence
        os.system("mv {} outputs/".format(name))

    """ Parse an orthogroup list to keep only one paralog sequence per species & per group
        (Keeps the 1st paralogous encoutered). No stringent, used in previous version
        Returns a list """
    def filterParalogsSoft(list_orthogroups, minspec):
        list_orthogroups_format = []
        j = 1

        for nofilter_group in list_orthogroups:
            new_group = []
            species = {}
            for loci in nofilter_group:
                species[loci.getHeader()[1:3]] = False
            for loci in nofilter_group:
                if not species[loci.getHeader()[1:3]]:
                    new_group.append(loci)
                    species[loci.getHeader()[1:3]] = True

            if len(new_group) >= minspec: # Drop too small orthogroups        
                list_orthogroups_format.append(new_group)
                writeOutputFile(new_group, j)
                j += 1
                
        return list_orthogroups_format

    """ Parse an orthogroup list to keep only the groups without any paralogs sequences
        Returns a list """
    def filterParalogsStringent(list_orthogroups, minspec):
        list_orthogroups_format = []
        j = 1

        for group in list_orthogroups:            
            species = {}
            for loci in group:
                species[loci.getHeader()[1:3]] = 0
            for loci in group:
                species[loci.getHeader()[1:3]] += 1            
            if len(set(species.values())) == 1 and 1 in set(species.values()) and len(group) >= minspec:                
                list_orthogroups_format.append(group)
                writeOutputFile(group, j)
                j += 1

        return list_orthogroups_format

    """ Builds a 2D array for a summary
        Returns a numpy 2D array """
    def countings(listOrthogroups, nb_rbh):

        def compute_nbspec(nb_rbh):
    
            def factorielle(x):
                n = 1
                s = 0
                while n <= x:
                    s += n
                    n += 1
                return s
            
            x = 2    
            nb_specs = 0
            while x*x - factorielle(x) < nb_rbh:
                x += 1 
            return x
        #listOrthogroups.sort().reverse()
        #nblines = len(listOrthogroups[0])
        nblines = 0    
        for group in listOrthogroups:
            if len(group) > nblines:
                nblines = len(group)
        matrix = np.array([[0]*compute_nbspec(nb_rbh)]*nblines)        

        for group in listOrthogroups:
            listSpecs = []
            for loci in group:
                if loci.getHeader()[1:3] not in listSpecs:
                    listSpecs.append(loci.getHeader()[1:3])
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
        return df

    # Function -------------------------------------------------------------------------------------------------
    list_orthogroups = []

    ## main loop : gather all pairwise couples who share loci ("My friend is your friend" : transitivity). They are then 
    #  written in fasta files.
    #  LAST EDIT : Galaxy can't handle properly collections with thousands of files. I applied some filters after this loop.
    #   in order to reduce the number of files returned. I could have modified the main loop instead, to reduce computation 
    #   time and memory consumption, but I have chosen to let it like this in the case Galaxy is improved in the future.
    #   (it will be easier to go back to the previous script).
    for ortho_pair1 in list_pairwises_allsp[0:-1]:
        if not checkIfTagged(ortho_pair1):
            orthogroup = ortho_pair1 # the orthogroup grows as we go throught the second loop

            # check for common locus between two groups
            for ortho_pair2 in list_pairwises_allsp[list_pairwises_allsp.index(ortho_pair1) + 1:]:
                if len(orthogroup.intersection(ortho_pair2)) != 0 and not checkIfTagged(ortho_pair2):
                    orthogroup.update(orthogroup | ortho_pair2)
                    tagGroup(ortho_pair2)

            # Check if subgroup is already computed
            if len(list_orthogroups) > 0: # every iterration except teh 1st one
                presence = False
                for group in list_orthogroups:
                    if len(group.intersection(orthogroup)) != 0:
                        group.update(group | orthogroup)
                        presence = True
                if not presence:
                    list_orthogroups.append(orthogroup)
            else:
                list_orthogroups.append(orthogroup) # necessary for the first iteration, when list_orthogroups is still empty
    
    # Options --------------------------------------------------------------------------------------------------

    """ nb : I could try to implement a more complex code which does in the same previous loop all the following lines, to avoid multiples parsing of
        the orthogroups list, but the code would become hardly readable. Since the whole program is relatively fast, I chosed code simplicity
        over code efficiency """

    # Print summary table with all paralogs
    if verbose :
        frame = countings(list_orthogroups, nb_rbh)
        df = asFrame(frame)
        print "\n    Summary before paralogous filtering : \n"
        print df.loc[df.ne(0).any(1),df.ne(0).any()], "\n" # Don't display columns and lines filled with 0    

    # Paralogs filtering and summary ----------------------------------------------------------------------------

    print "Filtering paralogous sequences and writing final orthogroups files ..."
    print "    (Dropping Orthogroups with less than {} species)".format(minspec)

    # writeOutputFile() is called in filterParalogsSoft()    
    list_orthogroups_format = filterParalogsStringent(list_orthogroups, minspec)   

    frame = countings(list_orthogroups_format, nb_rbh)

    df = asFrame(frame)
    print "\n    Summary after paralogous filtering : \n"
    print df.loc[df.ne(0).any(1),df.ne(0).any()]

    #return only the length of the list (at this point the program doesn't need more)
    return len(list_orthogroups_format)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("files", help="Input files separated by commas. Each file contains all the reciprocical best hits between a pair of species")
    parser.add_argument("minspec", help="Only keep Orthogroups with at least this number of species", type=int)    
    parser.add_argument("-v", "--verbose", action="store_true", help="A supplemental summary table of orthogroups before paralogs filtering will be returned")    
    args = parser.parse_args()

    print "*** pogs.py ***"
    print "\nBuilding of orthogroups based on pairs of genes obtained by pairwise comparisons between pairs of species."
    print "Genes are gathered in orthogroups based on the principle of transitivity between genes pairs."    

    os.system("mkdir outputs")
    
    infiles = args.files
    listPairwiseFiles = str.split(infiles, ",")
    print "\nParsing input files ..."
    list_Locus = getListPairwiseAll(listPairwiseFiles)
    print "Creating Orthogroups ..."    
    nb_orthogroups = makeOrthogroups(list_Locus, args.minspec, len(listPairwiseFiles), args.verbose)
    print "\n{} orthogroups have been infered from {} pairwise comparisons by RBH\n".format(nb_orthogroups, len(listPairwiseFiles))

if __name__ == "__main__":
    main()
