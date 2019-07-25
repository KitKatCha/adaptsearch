#!/usr/bin/env python
#coding: utf-8
#Author : Eric Fontanillas (2010) - Victor Mataigne (2018)

import pandas as pd
import argparse, os

def loop_on_elems(list_of_elems, path_in, path_out, sps_group_1, sps_group_2, colnames):

    # sub-routine

    def tableu(fileu, sps_group_1, sps_group_2):
        """ a

        Args :
            fileu : input file with counts of AAs / AAs types per orthogorup
            sps_group_1 : species for condition 1 (ex : hot water species)
            sps_group_2 : species for condition 2 (ex : cold water species)

        Returns :
            greater_dict :
            lower_dict :

        """
        df = pd.read_csv(fileu, sep=',', index_col=0, header=0)
        #print df
        # species = list(df) #columns names = species names

        # initialize counts
        greater_dict = {}
        lower_dict = {}

        for specie in sps_group_1+sps_group_2:
            greater_dict[specie] = 0
            lower_dict[specie] = 0

        #nb_trials = 0
        for (index, row) in df.iterrows():
            # min and max counts for each condition
            if not df.loc[index, sps_group_1+sps_group_2].isnull().values.any() :

                #nb_trials += 1

                max_cat1 = max(df.loc[index, sps_group_1]) # species in category 1 (ex : hots)
                min_cat1 = min(df.loc[index, sps_group_1])
                max_cat2 = max(df.loc[index, sps_group_2]) # species in category 2 (ex : colds)
                min_cat2 = min(df.loc[index, sps_group_2])

                for specie in sps_group_1:
                    if df.loc[index, specie] > max_cat2 :
                        greater_dict[specie] += 1
                    elif df.loc[index, specie] < min_cat2 :
                        lower_dict[specie] += 1

                for specie in sps_group_2:
                    if df.loc[index, specie] > max_cat1 :
                        greater_dict[specie] += 1
                    elif df.loc[index, specie] < min_cat1 :
                        lower_dict[specie] += 1

        return greater_dict, lower_dict#, nb_trials

    # Function ------------------------------------------------------

    for variable in list_of_elems:
        print 'Processing : {} ...'.format(variable)
        file_in = "{}/{}.csv".format(path_in, variable)
        file_out = open('{}/{}.csv'.format(path_out,variable), 'w')

        # Compute succeses and fails on each variable
        greater_dict, lower_dict = tableu(file_in, sps_group_1, sps_group_2)

        # totals and diffs
        diff_dict = {}
        total_dict = {}
        for key in greater_dict.keys():
            diff_dict[key] = greater_dict[key] - lower_dict[key]
            total_dict[key] = greater_dict[key] + lower_dict[key]
            #total_dict[key] = number_trials

        # results frame
        df = pd.DataFrame([greater_dict, lower_dict, diff_dict, total_dict])
        df = df.rename({0:'Greater',1:'Lower',2:'Difference',3:'Trial_Number'}) #, axis='index' if pandas 0.15
        df = df.rename(index=str, columns=colnames)

        df.to_csv("{}/{}.csv".format(path_out, variable), sep=",", encoding="utf-8")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sps_group_1", help="List of species separated by commas")
    parser.add_argument("sps_group_2", help="List of species separated by commas")
    parser.add_argument("format", choices=['nucleic', 'proteic'], help="input files format")
    args = parser.parse_args()

    # used only if format = nucleic
    LN =['A','C','T','G']
    Lratios = ['GC_percent', 'purine_percent', 'DIFF_GC', 'DIFF_AT', 'PLI_GC', 'PLI_AT', 'PLI_GC_1000', 'PLI_AT_1000']

    # used only if format = proteic
    LAA =['K','R','A','F','I','L','M','V','W','N','Q','S','T','H','Y','C','D','E','P','G']
    LV = ['IVYWREL','EK','ERK','DNQTSHA','QH','ratio_ERK_DNQTSHA','ratio_EK_QH','FYMINK','GARP',
          'ratio_GARP_FYMINK','AVLIM','FYW','AVLIMFYW','STNQ','RHK','DE','RHKDE','APGC','AC',
          'VLIM','ratio_AC_VLIM','ratio_APGC_VLIM']
    LS = ['total_residue_weight', 'total_residue_volume', 'total_partial_specific_volume', 'total_hydratation']

    # inputs and outputs paths
    if args.format == 'nucleic':
        input_path_elem = '02_tables_per_nucleotide'
        input_path_var = '02_tables_per_nuc_variable'
        out_path_elem = '03_tables_counts_signTest_nucleotides'
        out_path_var = '03_tables_counts_signTest_nuc_variables'
    elif args.format == 'proteic':
        input_path_elem = '02_tables_per_aa'
        input_path_var = '02_tables_per_aa_variable'
        out_path_elem = '03_tables_counts_signTest_aa'
        out_path_var = '03_tables_counts_signTest_aa_variables'

    os.mkdir(out_path_elem)
    os.mkdir(out_path_var)

    sps_group_1 = args.sps_group_1.split(',')
    sps_group_2 = args.sps_group_2.split(',')

    # Prepare colnames for final frames
    colnames = {}
    # for specie in sps_group_1:
    #     colnames[specie] = '{}_vs_condition_1'.format(specie)
    # for specie in sps_group_2:
    #     colnames[specie] = '{}_vs_condition_2'.format(specie)

    for specie in sps_group_1:
        colnames[specie] = '{}_VS_{}'.format(specie, args.sps_group_2.replace(',','_&_'))
    for specie in sps_group_2:
        colnames[specie] = '{}_VS_{}'.format(specie, args.sps_group_1.replace(',','_&_'))

    # Building tables
    if args.format == 'nucleic':
        loop_on_elems(LN, input_path_elem, out_path_elem, sps_group_1, sps_group_2, colnames)
        loop_on_elems(Lratios, input_path_var, out_path_var, sps_group_1, sps_group_2, colnames)
    elif args.format == 'proteic':
        loop_on_elems(LAA, input_path_elem, out_path_elem, sps_group_1, sps_group_2, colnames)
        loop_on_elems(LV, input_path_var, out_path_var, sps_group_1, sps_group_2, colnames)
        loop_on_elems(LS, input_path_var, out_path_var, sps_group_1, sps_group_2, colnames)

    # Final R script launching sign test
    print 'Processing : binomial sign tests ...'

    if args.format == 'nucleic':
        final_output_elem = '04_outputs_nucleotides'
        final_output_var = '04_outputs_nuc_variables'
    elif args.format == 'proteic':
        final_output_elem = '04_outputs_aa'
        final_output_var = '04_outputs_aa_variables'

    os.mkdir(final_output_elem)
    os.mkdir(final_output_var)
    os.system('Rscript S03b_sign_test_binomial.R --indir %s --outdir %s' %(out_path_elem, final_output_elem))
    os.system('Rscript S03b_sign_test_binomial.R --indir %s --outdir %s' %(out_path_var, final_output_var))

if __name__ == '__main__':
    main()
