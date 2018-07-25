#!/usr/bin/env python
# coding: utf-8
# Author : Victor Mataigne

from ete3 import EvolTree
from ete3 import Tree
#from ete3.treeview.layouts import evol_clean_layout #PyQt5
import argparse, csv, numpy
import pandas as pd
from functions_positive_selection import print_sites, models_types, frame_sites

# ./compCodeML_ete_oneFile.py locus_1_sp5.fasta tree.nwk M0,M1,M10,M11,M12,M13,M2,M3,M4,M5,M6,M7,M8a,M9

# Activate conda environment on nz cluster at Roscoff : source $CONDA3/activate ete3-3.1.1
    # but no pandas !
# Activate conda environment on vmataigne : source activate ete3

# TODOs:
#   - Make 'triangle' pandas dataFrames to avoid double NaN
#   - Finish the function to find positive selection
#   - Resolve dependancies/versions conflicts for visualization on GUI or in png.
#   - Make it callable in another script to run parallel computing on a set of input alignments
#   - For now, ML ratios comparisons are not done on models which have the same name. 
#     This is not a problem in the other script, 'load_precomputed_models', because the user can choose the model name, but in this script models are named according to
#     their built-in ete3 name

# This script performs codeML analysis on a alignment, using the tree file computed by Concatphyl (super-alignment) or, if available, on the tree for this alignment only.
# It computes all the evolutionary models specified in the models_list command-line argument

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('Alignment', help='Alignment fasta file')
    parser.add_argument('Tree', help='Phylogenetic tree in newick format')
    parser.add_argument('models_list', help='List of evolutionary models to run')
    parser.add_argument('-v', '--verbose', action="store_true", help='Display infos on conserved sites positions')
    args = parser.parse_args()
    
    # Models list
    inputModels = args.models_list
    models_list = inputModels.split(',')

    # Set tree
    inputTree = open(args.Tree, 'r')
    nwTree = inputTree.readline()
    inputTree.close()
    tree = EvolTree(nwTree)

    # Set alignment
    inputAlig = open(args.Alignment, 'r')
    alig = inputAlig.read()
    inputAlig.close()

    # Run selected models on each alignment
    print '\nPart I - Running all models\n---------------------------\n'

    if alig.count('>') == len(tree.get_leaves()):    
        tree.link_to_alignment(alig)
        for model in models_list:
            #if model in ete_available_models:
            print '        Running model {} ...'.format(model)
            tree.run_model(model)
            #else :
                #print '    Error with {} : unknown model'.format(model)
    else :
        print '    Error in {} : number of species not equal to number of leaves'.format(args.Alignment)
        raise SystemExit

    # Display models on screen / in log file
    print ""
    for model in tree._models:
        print tree.get_evol_model(model)

    # Testing models against each other
    print '\nPart II - Testing models\n------------------------\n'           

    # argument alt: model with higher number of parameters (np)
    # argument null: model with lower number of parameters (np)
    
    table = numpy.empty((len(tree._models),len(tree._models)))
    table[:] = numpy.NaN
    df = pd.DataFrame(data=table, columns=models_list, index=models_list)
    df = df.rename_axis('Alt. model').rename_axis('Null model', axis='columns')

    for alt_model in tree._models:
        alt = tree.get_evol_model(alt_model)
        for null_model in tree._models:
            null = tree.get_evol_model(null_model)
            if alt_model != null_model and alt.np > null.np and null.lnL - alt.lnL < 0:
                pval = tree.get_most_likely(alt_model, null_model)
                df.at[alt_model, null_model] = pval
            elif null.lnL - alt.lnL > 0 and alt.np > null.np: 
                df.at[alt_model, null_model] = 1

    
    print "Summary table :\n"
    print "    The test hypothesis is that the 'Alternative model' fits more than the null model.\n"
    print "    - NaN values mean that hypothesis testing was impossible because the null model has more paramaters than the alternative model:"
    print "    - 1 values mean that hypothesis testing was impossible because the log-likelihood of the alternative model is smaller than null's"
    print "    pvalue < 0.05 means that the alternative model wins.\n"
    print df, '\n'

    #tree.show(histfaces=['M2'])

    #tree.render(file_name=args.Alignment.split('.')[0], layout=evol_clean_layout, histfaces=['M1'])
    print "Part III - Detecting selection\n------------------------------\n"
    mtypes = models_types()
    for model in tree._models:
        print_sites(model, tree, mtypes, args.verbose) # print on screen/log file
        frame_sites(model, tree, mtypes, args.verbose) # write in csv table

if __name__ == '__main__':
    main()

 # The config file looks like this :

 #      seqfile = fileName * sequence data file name
 #      outfile = nameResults * main result file name
 #     treefile = treeName * tree structure file name
 #        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
 #      verbose = 0  * 1: detailed output, 0: concise output
 #      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
 #                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
 #      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
 #    CodonFreq = 1  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
 #        clock = 0   * 0:no clock, 1:clock; 2:local clock
 #       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
 #                   * 7:AAClasses
 #   aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)
 #                  * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
 #        model = 0
 #                   * models for codons:
 #                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
 #                   * models for AAs or codon-translated AAs:
 #                       * 0:poisson, 1:proportional,2:Empirical,3:Empirical+F
 #                       * 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)
 #      NSsites = 1  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
 #                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
 #                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
 #                   * 13:3normal>0
 #        icode = 0  * 0:universal code; 1:mammalian mt; 2-11:see below
 #        Mgene = 0  * 0:rates, 1:separate;
 #    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
 #        kappa = 2.0  * initial or fixed kappa
 #    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
 #        omega = 0.2 * initial or fixed omega, for codons or codon-based AAs
 #    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
 #        alpha = 0.0 * initial or fixed alpha, 0:infinity (constant rate)
 #       Malpha = 0  * 1: different alphas for genes, 0 : one alpha
 #        ncatG = 3  * # of categories in dG of NSsites models
 #      fix_rho = 1  * 0: estimate rho; 1: fix it at rho
 #          rho = 0.0 * initial or fixed rho,   0:no correlation
 #        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 # RateAncestor = 0  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
 #   Small_Diff = 5e-07
 #    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?
 #  fix_blength = 0   * 0: ignore, -1: random, 1: initial, 2: fixed
 #       method = 0   * 0: simultaneous; 1: one branch at a time

# Available models

# Name   | Description                 | Model kind
# ------------------------------------------------
# M1     | relaxation                  | site
# M10    | beta and gamma + 1          | site
# M11    | beta and normal > 1         | site
# M12    | 0 and 2 normal > 2          | site
# M13    | 3 normal > 0                | site
# M2     | positive-selection          | site
# M3     | discrete                    | site
# M4     | frequencies                 | site
# M5     | gamma                       | site
# M6     | 2 gamma                     | site
# M7     | relaxation                  | site
# M8     | positive-selection          | site
# M8a    | relaxation                  | site
# M9     | beta and gamma              | site
# SLR    | positive/negative selection | site
# M0     | negative-selection          | null
# fb_anc | free-ratios                 | branch_ancestor
# bsA    | positive-selection          | branch-site
# bsA1   | relaxation                  | branch-site
# bsB    | positive-selection          | branch-site
# bsC    | different-ratios            | branch-site
# bsD    | different-ratios            | branch-site
# b_free | positive-selection          | branch
# b_neut | relaxation                  | branch
# fb     | free-ratios                 | branch
# XX     | User defined                | Unknown

# Usual comparisons are :

# ============ ======= ===========================================
#  Alternative  Null    Test
# ============ ======= ===========================================
#   M2          M1      PS on sites (M2 prone to miss some sites)
#                       (Yang 2000).
#   M3          M0      test of variability among sites
#   M8          M7      PS on sites
#                       (Yang 2000)
#   M8          M8a     RX on sites?? think so....
#   bsA         bsA1    PS on sites on specific branch
#                       (Zhang 2005)
#   bsA         M1      RX on sites on specific branch
#                       (Zhang 2005)
#   bsC         M1      different omegas on clades branches sites
#                       ref: Yang Nielsen 2002
#   bsD         M3      different omegas on clades branches sites
#                       (Yang Nielsen 2002, Bielawski 2004)
#   b_free      b_neut  foreground branch not neutral (w != 1)
#                        - RX if P<0.05 (means that w on frg=1)
#                        - PS if P>0.05 and wfrg>1
#                        - CN if P>0.05 and wfrg>1
#                        (Yang Nielsen 2002)
#   b_free      M0      different ratio on branches
#                       (Yang Nielsen 2002)
# ============ ======= ===========================================
# **Note that M1 and M2 models are making reference to the new versions
# of these models, with continuous omega rates (namely M1a and M2a in the
# PAML user guide).**

# **Alternative must have a greater number of parameters than Null