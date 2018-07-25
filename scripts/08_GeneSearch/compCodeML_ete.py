#!/usr/bin/env python
# coding: utf-8
# Author : Victor Mataigne

from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('Alignments', help='Alignment fasta files')
    parser.add_argument('Tree', help='Phylogenetic tree in newick format')    
    #parser.add_argument('configFile', help='codeml.ctl file')
    parser.add_argument('models_list', help='List of evolutionary models to run')
    args = parser.parse_args()

    # Alignment files
    infiles = args.Alignments
    inputs = infiles.split(',') # list of input files (if one file : list of 1)

    # Models list
    inputModels = args.models_list
    models_list = inputModels.split(',')

    # Set tree
    inputTree = open(args.Tree, 'r')
    nwTree = inputTree.readline()
    inputTree.close()
    

    # Store every models for every alignment
    models_per_alig = {}    

    # Run selected models on each alignment
    print 'Part I - Running all models'
    for file in inputs :

        # Cannot create the tree outside of the loop because it keeps all links on all alignments, so same models are computed several times
        tree = EvolTree(nwTree)

        # Set alignment
        inputAlig = open(file, 'r')
        alig = inputAlig.read()
        inputAlig.close()

        if alig.count('>') == len(tree.get_leaves()):
        
            tree.link_to_alignment(alig)

            print '    Running models on {} ...'.format(file)

            for model in models_list:
                #if model in ete_available_models:
                print '        Running model {} ...'.format(model)
                tree.run_model(model)
                #else :
                    #print '    Error with {} : unknown model'.format(model)

            models_per_alig[file] = tree

        else :
            print '    Error in {} : number of species not equal to number of leaves'.format(file)

    # Display models on screen / in log file
    for key in models_per_alig.keys():
        print '\nAlignment : {}'.format(key)
        print '    Models :'
        for model in models_per_alig[key]._models:
            print models_per_alig[key].get_evol_model(model)

    # Testing models against each other
    print 'Part II - Testing models'
    
    for key in models_per_alig.keys():        
        print "\n    {}\n".format(key)
        i = 0
        table = np.zeros((len(models_per_alig[key]._models),len(models_per_alig[key]._models)))
        for alt_model in models_per_alig[key]._models:            
            j = 0
            for null_model in models_per_alig[key]._models:
                pval = tree.get_most_likely(alt_model, null_model)
                table[i][j] = pval
                j += 1
            i += 1

        print table

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

# Hypothesis testing

# | Alternative | Null   | Test
# --------------------------------------------------------------------------------------------------------------------
# |   M2        | M1     |    PS on sites (M2 prone to miss some sites) (Ziheng Yang (2000))
# |   M3        | M0     |    test of variability among sites
# |   M8        | M7     |    PS on sites (Ziheng Yang (2000))
# |   M8        | M8a    |    RX on sites
# |   bsA       | bsA1   |    PS on sites on specific branch (Jianzhi Zhang (2005))
# |   bsA       | M1     |    RX on sites on specific branch (Jianzhi Zhang (2005))
# |   bsC       | M1     |    different omegas on clades branches sites (Ziheng Yang (2002))
# |   bsD       | M3     |    different omegas on clades branches sites (Ziheng Yang (2002), Joseph P. Bielawski (2004))
# |   b_free    | b_neut |    foreground branch not neutral (w != 1)
# |             |        |    - RX if P<0.05 (means that w on frg=1)
# |             |        |    - PS if P>0.05 and wfrg>1
# |             |        |    - CN if P>0.05 and wfrg>1
# |             |        |    (Ziheng Yang (2002))
# |   b_free    | M0     |    different ratio on branches (Ziheng Yang (2002))
