#!/usr/bin/env python
#coding: utf-8

import argparse, os, re
from functions import dico, write_output, fill_with_NaN

""" Functions for nucleic format """

def all_nuc_counts(seq):
    """ Counts the number of nucleotides in a sequence

    Args :
        - seq (String) : a nucleic sequence

    Return:
        - nuc_counts (dict) : counts on nucleotides (key/value : nucleotide/count)
    """
    nuc_counts = {}
    seqU = seq.upper()
    LN =['A','C','T','G']

    for base in LN:
        nuc_counts[base] = seqU.count(base)

    return nuc_counts

def ratios(nuc_counts):
    """ Compute GC and purine ratios from the counts of nucleotides in a sequence

    Args :
        - nuc_counts (dict) : dictionary of nucleotides counts (output of all_nuc_counts())

    Return :
        - ratios (dict) : dictionary with ratios
    """

    # Search for compositional bias in genome as marker of thermal adaptation
    ratios = {}
    ratios['GC_percent'] = float(nuc_counts['C'] + nuc_counts['G'])/(sum(nuc_counts.values())*100)
    ratios['purine_percent'] = float(nuc_counts['A'] + nuc_counts['G'])/(sum(nuc_counts.values())*100)

    ratios['DIFF_GC'] = nuc_counts['G'] - nuc_counts['C']
    ratios['DIFF_AT'] = nuc_counts['A'] - nuc_counts['T']

    # Per bp
    ratios['PLI_GC'] = float(ratios['DIFF_GC'])/sum(nuc_counts.values())
    ratios['PLI_AT'] = float(ratios['DIFF_AT'])/sum(nuc_counts.values())

    # Per 1000 bp
    ratios['PLI_GC_1000'] = ratios['PLI_GC']*1000
    ratios['PLI_AT_1000'] = ratios['PLI_AT']*1000

    return ratios

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("species_list", help="List of species separated by commas")
    args = parser.parse_args()

    LN =['A','C','T','G']
    Lratios = ['GC_percent', 'purine_percent', 'DIFF_GC', 'DIFF_AT', 'PLI_GC', 'PLI_AT', 'PLI_GC_1000', 'PLI_AT_1000']

    list_inputs = []

    path_inputs = '01_input_files'
    list_inputs = os.listdir(path_inputs)

    lsp = args.species_list.split(',')
    lsp = sorted(lsp)
    flsp = ''
    for el in lsp:
        flsp += el+','

    path_outputs_1 = '02_tables_per_nucleotide'
    path_outputs_2 = '02_tables_per_nuc_variable'
    os.mkdir(path_outputs_1)
    os.mkdir(path_outputs_2)

    dict_for_files_nuc = {}
    dict_for_file_var = {}

    # All counts
    for file in list_inputs:
        # iterate over input files
        sequences = dico(file, path_inputs, lsp)

        # TEMPORARY CORRECTION FOR SEQUENCES CONTAINING ONLY INDELS
        # It appears than CDS_Search can bug sometimes and return an alignement where a species' sequence is made of indels only
        # This causes a crash here (in the ratios function). The correction skip the whole file.
        # When CDS_Search is corrected, lines with 'skip' can be removed
        skip = False
        for key in sequences.keys():
            if all(x == '-' for x in sequences[key]):
                skip = True

        if not skip:

            nuc_counts_per_seq = {}
            nuc_var_per_seq = {}

            for key in sequences.keys():
                # iterate over sequences in the file
                nuc_counts_per_seq[key] = all_nuc_counts(sequences[key])
                nuc_var_per_seq[key] = ratios(nuc_counts_per_seq[key])

            # Add NaN for missing species
            for key in set(lsp).difference(set(sequences.keys())):
                nuc_counts_per_seq[key] = fill_with_NaN(LN)
                nuc_var_per_seq[key] = fill_with_NaN(Lratios)

            # Add computations to final dicts
            dict_for_files_nuc[file] = nuc_counts_per_seq
            dict_for_file_var[file] = nuc_var_per_seq


    # Try with pandas ?
    write_output(LN, flsp, path_outputs_1, dict_for_files_nuc) # one file per nuc
    write_output(Lratios, flsp, path_outputs_2, dict_for_file_var) # one file per nuc_variable

if __name__ == '__main__':
    main()
