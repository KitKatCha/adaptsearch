#!/usr/bin/env python
#coding: utf-8
#Author : Eric Fontanillas (2010) - Victor Mataigne (2018)

# TODO :
    # - Deal with missing data : do not do the sign test if missing species in a group
    # - Find a way to avoid the list_species argument

import argparse, os, re
from functions import dico, write_output, fill_with_NaN

def aa_properties(amino_acids_properties_file):
    """ Read the file 'amino_acids_properties' and stores its content

    Args :
        amino_acids_properties_file (String) : the file

    Return :
        aa_properties (dict) : key/amino-acid - value/list of properties
    """

    dict_aa_properties={}
    with open(amino_acids_properties_file, 'r') as f:
        f.readline() #jump headers
        for line in f.readlines():
            S1 = line.split(",")
            aa_name = S1[1]
            S2 = aa_name.split("/")
            aa_code = S2[1][:-1]

            frequencies = S1[2][:-1]
            residue_weight = S1[5]
            residue_volume = S1[6]
            partial_specific_volume = S1[7]
            hydration = S1[8]

            dict_aa_properties[aa_code] = [frequencies, residue_weight, residue_volume, partial_specific_volume, hydration]

    return(dict_aa_properties)

""" Functions for proteic format """

def all_aa_counts(seq):
    """ Count the occurrences of all amino-acids in a sequence

    Args:
        seq (String) : a proteic sequence

    Returns: a dictionary with amino-acids counts
    """

    aa_counts = {}
    seqU = seq.upper()
    LAA =['K','R','A','F','I','L','M','V','W','N','Q','S','T','H','Y','C','D','E','P','G']

    for aa in LAA:
        aa_counts[aa] = seqU.count(aa)

    return aa_counts

def all_aa_props(seq_counts):
    """ Converts a dictionnary of counts into a dictionnary of proportions

    Args:
        seq_counts (dict) : dictionnary computed by the function all_aa_counts()

    Returns: a dictionary with counts replaced by proportions
    """

    aa_props = {}
    for key in seq_counts.keys():
        aa_props[key] = float(seq_counts[key]) / sum(seq_counts.values())
    return aa_props

def aa_variables_counts_and_props(aa_counts):
    """ Computes several thermostability indices (summed occurrences of some AAs, and then various ratios)

    Args:
        aa_counts (dict) : dictionnary computed by the function all_aa_counts()

    Returns:
        aa_variables_counts : a dictionary with indices values
        aa_variables_props : a dictionary with indices proportions (ratios excluded)
    """

    # Hyperthermophile Prokaryotes criterias
    #   IVYWREL : positivelly correlated with otpimal growth
    #   ERK : (i.e. ) => positivelly correlated with optimal growth temperature
    #     ERK/DNQTSHA (or DNQTSH ??)
    #     EK/QH

    # Mutationnal bias hypothesis => AT rich: favor FYMINK // GC rich: favor GARP
    #      The mutational bias model predict a linear relationship between GARP vs FYMINK
    #      ==> so if outliers to that, it means that the excess of GARP or FYMINK are not explained by the mutationnal bias model but by something else (selection ?)


    # Hydophobicity hypothesis [should INCREASE with thermal adaptation]
    #   AL
    #   Only non-aromatic : AVLIM
    #   Only aromatic : FYW

    # Charged hypothesis => positivelly correlated with optimal growth temperature
    #   All charged : RHKDE
    #   Only positive : RHK
    #   Only negative : DE

    # Neutral polar hypothesis  [should DECREASE with thermal adaptation]
    #   STNQ

    # Fontanillas' criteria
    #   PAYRE
    #   MVGDS
    #   PAYRE/MVGDS

    # Jollivet's criteria
    #   AC
    #   MVGDS
    #   AC/MVGDS

    aa_variables_counts = {}
    aa_variables_props = {}
    len_seq = sum(aa_counts.values()) # length of the sequence

    # counts of variables
    aa_variables_counts['AC'] = aa_counts['A'] + aa_counts['C']
    aa_variables_counts['APGC'] = aa_counts['A'] + aa_counts['P'] + aa_counts['G'] + aa_counts['C']
    aa_variables_counts['AVLIM'] = aa_counts['A'] + aa_counts['V'] + aa_counts['L'] + aa_counts['I'] + aa_counts['M']
    aa_variables_counts['AVLIMFYW'] = aa_variables_counts['AVLIM'] + aa_counts['F'] + aa_counts['Y'] + aa_counts['W']
    aa_variables_counts['DE'] = aa_counts['D'] + aa_counts['E']
    aa_variables_counts['DNQTSHA'] = aa_counts['D'] + aa_counts['N'] + aa_counts['Q'] + aa_counts['T'] + aa_counts['S'] + aa_counts['H'] + aa_counts['A']
    aa_variables_counts['EK'] = aa_counts['E'] + aa_counts['K']
    aa_variables_counts['ERK'] = aa_counts['E'] + aa_counts['K'] + aa_counts['K']
    aa_variables_counts['FYMINK'] = aa_counts['F'] + aa_counts['Y'] + aa_counts['M'] + aa_counts['I'] + aa_counts['N'] + aa_counts['K']
    aa_variables_counts['FYW'] = aa_counts['F'] + aa_counts['Y'] + aa_counts['W']
    aa_variables_counts['GARP'] = aa_counts['G'] + aa_counts['A'] + aa_counts['R'] + aa_counts['P']
    aa_variables_counts['IVYWREL'] = aa_counts['I'] + aa_counts['V'] + aa_counts['Y'] + aa_counts['W'] + aa_counts['R'] + aa_counts['E'] + aa_counts['L']
    aa_variables_counts['QH'] = aa_counts['Q'] + aa_counts['H']
    aa_variables_counts['RHK'] = aa_counts['R'] + aa_counts['H'] + aa_counts['K']
    aa_variables_counts['RHKDE'] = aa_counts['R'] + aa_counts['H'] + aa_counts['K'] + aa_counts['D'] + aa_counts['E']
    aa_variables_counts['STNQ'] = aa_counts['S'] + aa_counts['T'] + aa_counts['N'] + aa_counts['Q']
    aa_variables_counts['VLIM'] = aa_counts['V'] + aa_counts['L'] + aa_counts['I'] + aa_counts['M']
    aa_variables_counts['PAYRE'] = aa_counts['P'] + aa_counts['A'] + aa_counts['Y'] + aa_counts['R'] + aa_counts['E']
    aa_variables_counts['MVGDS'] = aa_counts['M'] + aa_counts['V'] + aa_counts['G'] + aa_counts['D'] + aa_counts['S']

    # compute proportions
    for key in aa_variables_counts.keys():
        aa_variables_props[key] = float(aa_variables_counts[key]) / float(len_seq)

    if aa_variables_counts['DNQTSHA'] != 0:
        ratio_ERK_DNQTSHA = float(aa_variables_counts['ERK'])/float(aa_variables_counts['DNQTSHA'])
    else :
        ratio_ERK_DNQTSHA = -1

    if aa_variables_counts['QH'] != 0:
        ratio_EK_QH = float(aa_variables_counts['EK'])/float(aa_variables_counts['QH'])
    else :
        ratio_EK_QH = -1

    if aa_variables_counts['FYMINK'] != 0:
        ratio_GARP_FYMINK = float(aa_variables_counts['EK'])/float(aa_variables_counts['FYMINK'])
    else :
        ratio_GARP_FYMINK = -1

    if aa_variables_counts['VLIM'] != 0:
        ratio_AC_VLIM = float(aa_variables_counts['AC'])/float(aa_variables_counts['VLIM'])
        ratio_APGC_VLIM =  float(aa_variables_counts['APGC'])/float(aa_variables_counts['VLIM'])
    else :
        ratio_AC_VLIM = -1
        ratio_APGC_VLIM = -1

    if aa_variables_counts['MVGDS'] != 0:
        ratio_PAYRE_MVGDS = float(aa_variables_counts['PAYRE'])/float(aa_variables_counts['MVGDS'])
    else :
        ratio_PAYRE_MVGDS = -1

    if aa_variables_counts['MVGDS'] != 0:
        ratio_AC_MVGDS = float(aa_variables_counts['AC'])/float(aa_variables_counts['MVGDS'])
    else :
        ratio_AC_MVGDS = -1

    aa_variables_counts['ratio_ERK_DNQTSHA'] = ratio_ERK_DNQTSHA
    aa_variables_counts['ratio_EK_QH'] = ratio_EK_QH
    aa_variables_counts['ratio_GARP_FYMINK'] = ratio_GARP_FYMINK
    aa_variables_counts['ratio_AC_VLIM'] = ratio_AC_VLIM
    aa_variables_counts['ratio_APGC_VLIM'] = ratio_APGC_VLIM
    aa_variables_counts['ratio_PAYRE_MVGDS'] = ratio_PAYRE_MVGDS
    aa_variables_counts['ratio_AC_MVGDS'] = ratio_AC_MVGDS

    return aa_variables_counts, aa_variables_props

def sequence_properties_from_aa_properties(aa_counts, aa_properties):
    """ Computes a sequence properties (based on an external data file)

    Args:
        - aa_counts (dict) : counts of amino-acids in the sequence
        - aa_properties (dict) : key/amino-acid - value/list of properties extract from the external data file

    Returns:
        - seq_props (dict) : values of the sequence properties
    """

    LS = ['total_residue_weight', 'total_residue_volume', 'total_partial_specific_volume', 'total_hydratation']
    seq_props = {}

    for i in range(1,5):
        seq_props[LS[i-1]] = 0
        for key in aa_counts.keys():
            seq_props[LS[i-1]] += aa_counts[key] * float(aa_properties[key][i])

    return seq_props

""" Main """

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("species_list", help="List of species separated by commas")
    parser.add_argument("aa_properties", help="File with all amino-acids properties")
    args = parser.parse_args()

    LAA = ['K','R','A','F','I','L','M','V','W','N','Q','S','T','H','Y','C','D','E','P','G']
    LV = ['IVYWREL','EK','ERK','DNQTSHA','QH','ratio_ERK_DNQTSHA','ratio_EK_QH','FYMINK','GARP',
          'ratio_GARP_FYMINK','AVLIM','FYW','AVLIMFYW','STNQ','RHK','DE','RHKDE','APGC','AC',
          'VLIM','ratio_AC_VLIM','ratio_APGC_VLIM']
    LS = ['total_residue_weight', 'total_residue_volume', 'total_partial_specific_volume', 'total_hydratation']

    list_inputs = []

    path_inputs = '01_input_files'
    list_inputs = os.listdir(path_inputs)

    lsp = args.species_list.split(',')
    lsp = sorted(lsp)
    flsp = ''
    for el in lsp:
        flsp += el+','

    path_outputs_1 = '02_tables_per_aa'
    path_outputs_2 = '02_tables_per_aa_variable'
    os.mkdir(path_outputs_1)
    os.mkdir(path_outputs_2)

    # Init empty dicts for results
    dict_for_files_aa_counts = {} # counts
    dict_for_files_aa_props = {} # proportions
    dict_for_files_variables_counts = {}
    dict_for_files_variables_props = {}
    dict_for_files_seq_properties = {}

    aa_properties_file = aa_properties(args.aa_properties) # read the aa_properties file

    # All counts and props
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

            aa_counts_per_seq = {}
            aa_props_per_seq = {}
            aa_variables_counts_per_seq = {}
            aa_variables_props_per_seq = {}
            seq_properties = {}

            for key in sequences.keys():
                # iterate over sequences in the file
                aa_counts_per_seq[key] = all_aa_counts(sequences[key])
                aa_props_per_seq[key] = all_aa_props(aa_counts_per_seq[key])
                aa_variables_counts_per_seq[key], aa_variables_props_per_seq[key] = aa_variables_counts_and_props(aa_counts_per_seq[key])
                seq_properties[key] = sequence_properties_from_aa_properties(aa_counts_per_seq[key], aa_properties_file)

            # Add NaN for missing species
            for key in set(lsp).difference(set(sequences.keys())):
                aa_counts_per_seq[key] = fill_with_NaN(LAA)
                aa_props_per_seq[key] = fill_with_NaN(LAA)
                aa_variables_counts_per_seq[key] = fill_with_NaN(LV)
                seq_properties[key] = fill_with_NaN(LS)

            # Add computations to final dicts
            dict_for_files_aa_counts[file] = aa_counts_per_seq
            dict_for_files_aa_props[file] = aa_props_per_seq
            dict_for_files_variables_counts[file] = aa_variables_counts_per_seq
            dict_for_files_variables_props[file] = aa_variables_props_per_seq
            dict_for_files_seq_properties[file] = seq_properties

    # Try with pandas ?
    write_output(LAA, flsp, path_outputs_1, dict_for_files_aa_counts) # one file per AA
    write_output(LV, flsp, path_outputs_2, dict_for_files_variables_counts) # one file per aa_variable
    write_output(LS, flsp, path_outputs_2, dict_for_files_seq_properties) #one file per seq properties

    # { file_name1 : { seq1: {'A' : 0, 'C':0, 'E':0, ...},
    #                  seq2: {'A' : 0, 'C':0, 'E':0, ...}, ...
    #                },
    # { file_name2 : { seq1: {'A' : 0, 'C':0, 'E':0, ...},
    #                  seq2: {'A' : 0, 'C':0, 'E':0, ...},
    #                },
    # ... }

    # { file_name1 : {seq1 : {'IVYWREL' : 0, 'EK': 0, 'ERK': 0, ...},
    #                 seq2 : {'IVYWREL' : 0, 'EK': 0, 'ERK': 0, ...}, ...
    #                },
    #   file_name2 : {seq1 : {'IVYWREL' : 0, 'EK': 0, 'ERK': 0, ...},
    #                 seq2 : {'IVYWREL' : 0, 'EK': 0, 'ERK': 0, ...},
    #                },
    #   ... }

    # { file_name1 : {'IVYWREL' : 0, 'EK': 0, 'ERK': 0, ...},
    #   file_name2 : {'IVYWREL' : 0, 'EK': 0, 'ERK': 0, ...},
    #   ...
    # }

if __name__ == '__main__':
    main()
