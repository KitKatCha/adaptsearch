#!/usr/bin/env python
# coding: utf8
# Author : Victor Mataigne

# ete3 evol -t RAxML_bestTree.nw --alg concat.fasta.fasta -o results1/ --models fb M2 SLR --cpu 3

# WARNING : 
#   - module ete3 : the one from conda envt named 'ete3' (conda2activate | source activate ete3)
#   - codeml binary : the one on nz
# -> work on nz while ete3 envt is activated

# Activate conda environment on nz cluster at Roscoff : source $CONDA3/activate ete3-3.1.1
    # but no pandas !
# Activate conda environment on vmataigne : source activate ete3

# TODOs:
#   - Make a 'triangle' pandas dataFrame to avoid double NaN
#   - Finish the function to find positive selection
#   - Resolve dependancies/versions conflicts for visualization on GUI or in png.

# http://etetoolkit.org/docs/latest/tutorial/tutorial_adaptation.html

from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
import argparse, numpy
import pandas as pd
from functions_positive_selection import details_most_likely_on_sites, frame_site

# --- main function is for reading codeML outputs and performing ML ratios tests.

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('Tree', help='Phylogenetic tree in newick format')
    parser.add_argument('Alignment', help='Alignment fasta files')
    parser.add_argument('models_list', help='File with list of pre-computed models. Each line must be as follows : <MODEL_NAME> <MODEL_PATH>')
    args = parser.parse_args()

    print "\n--- Testing Evolutionary Hypothesis with codeML and ETE ---\n"

    inputTree = open(args.Tree, 'r')
    nwTree = inputTree.readline()
    inputTree.close()
    inputAlig = open(args.Alignment, 'r')
    alig = inputAlig.read()
    inputAlig.close()    

    tree = EvolTree(nwTree)
    tree.link_to_alignment(alig)

    # Add pre-computed models 
    computed_models_names = []
    computed_models_paths = []
    models_names_display = ""

    with open(args.models_list, 'r') as read_models:
        for line in read_models.readlines():
            m_prop = line.split(' ')
            models_names_display += m_prop[0]+ ' '
            computed_models_names.append(m_prop[0])
            computed_models_paths.append(m_prop[1].strip('\n'))
            tree.link_to_evol_model(m_prop[1].strip('\n'), m_prop[0])

    # Run another model
    #tree.run_model('M0')

    print "Loaded models : %s\n" %models_names_display
    print "Models details :\n"
    for model in tree._models:
        print tree.get_evol_model(model)

    # Testing models against each other
    print 'Testing models :\n'

    # argument alt: model with higher number of parameters (np)
    # argument null: model with lower number of parameters (np)
    
    table = numpy.empty((len(tree._models),len(tree._models)))
    table[:] = numpy.NaN
    df = pd.DataFrame(data=table, columns=computed_models_names, index=computed_models_names)
    df = df.rename_axis('Alt. model').rename_axis('Null model', axis='columns')
    
    for alt_model in tree._models:
        alt = tree.get_evol_model(alt_model)
        for null_model in tree._models:
            null = tree.get_evol_model(null_model)
            # Checks if comparison is possible
            # the alt.name != null.name maybe problematic in the other script (in case of the comparison of two models of same type but with different parameters)
            if alt.name != null.name and alt.np > null.np and null.lnL - alt.lnL < 0:
                pval = tree.get_most_likely(alt_model, null_model)
                df.at[alt_model, null_model] = pval
    
    print "\nSummary table :\n"
    print "    The test hypothesis is that the 'Alternative model' fits more than the null model.\n"
    print "    NaN values mean that hypothesis testing was impossible either because :"
    print "        - likelihood of the alternative model is smaller than or equal to null's"
    print "        - null model has more paramaters than alternative's"
    print "        - alt and null models are the same (matrix diagonal top-left/bottom-right)"
    print "    pvalue < 0.05 means that the alternative model wins.\n"
    print df
    #print df.where(numpy.tril(numpy.ones(df.shape)).astype(numpy.bool))

    # Layout

    #tree.show(histfaces=['M2']) # 'QGraphicsSimpleTextItem' object has no attribute 'rotate'
    #tree.show(layout=evol_clean_layout)
    #tree.render(file_name='tree.png', layout=evol_clean_layout)
    #tree.render(file_name=args.Alignment.split('.')[0], layout=evol_clean_layout, histfaces=['M1'])

    # classical_comparisons = {'M2' : 'M1',
    #                 'M3' : 'M0',
    #                 'M8' : 'M7',
    #                 'M8' : 'M8a',
    #                 'bsA' : 'bsA1',
    #                 'bsA': 'M1',
    #                 'bsC' : 'M1',
    #                 'bsD' : 'M3',
    #                 'b_free' : 'b_neut',
    #                 'b_free' : 'M0'}

    # for model in classical_comparisons.keys():
    #     if model in tree._models and classical_comparisons[model] in tree._models:
    #         pval = tree.get_most_likely(model, classical_comparisons[model])
    #         # tests on sites
    #         if model.startswith('M') :                
    #             details_most_likely_on_sites(pval, model, classical_comparisons[model], tree)
    #         # tests on branch-sites
    #         elif model.startswith('bs') :
    #             details_most_likely_on_bs(pval, model, classical_comparisons[model], tree)
    #         # tests on branches
    #         elif model.startswith('b_'):
    #             details_most_likely_on_branches(pval, model, classical_comparisons[model], tree)

if __name__ == "__main__":
    main()

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