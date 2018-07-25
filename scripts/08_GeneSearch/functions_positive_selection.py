#!/usr/bin/env python
# coding: utf8
# Author : Victor Mataigne

from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
import argparse, numpy, csv
import pandas as pd

# -- Unfinished functions for finding details on positive selection

def print_sites(model, tree, models_types, verbose):

    def print_selection(Analysis, model, tree, verbose):
        print  '    - {} Analysis -\n'.format(Analysis)

        plist = {
            True: ['p0','p1','p2'],
            False: ['p2']
        }[verbose]

        for p in plist :
            if p in tree.get_evol_model(model).sites[Analysis]:
                print "\ndetails for probability {}\n".format(p)
                for s in range(len(tree.get_evol_model(model).sites[Analysis]['aa'])):
                    if tree.get_evol_model(model).sites[Analysis][p][s] > 0.95 and p == 'p2':
                        print 'positively selected site %s at position: %s, with probability: %s' %((tree.get_evol_model(model).sites[Analysis]['aa'][s], s+1, tree.get_evol_model(model).sites[Analysis][p][s]))
                    if tree.get_evol_model(model).sites[Analysis][p][s] > 0.95 and p == 'p1':
                        print 'Neutral site %s at position: %s, with probability: %s' %((tree.get_evol_model(model).sites[Analysis]['aa'][s], s+1, tree.get_evol_model(model).sites[Analysis][p][s]))
                    if tree.get_evol_model(model).sites[Analysis][p][s] > 0.95 and p == 'p0':
                        print 'Conserved site %s at position: %s, with probability: %s' %((tree.get_evol_model(model).sites[Analysis]['aa'][s], s+1, tree.get_evol_model(model).sites[Analysis][p][s]))

    if model in models_types["site"] :
        print "Model {} - List of positively selected sites :\n".format(model)

        print " For each of BEB and NEB analysis, the probability of belonging from one category of site is summarized by ‘p0’, ‘p1’ and" 
        print " ‘p2’ in the case of M2 model that have only 3 class of sites. Only positive-selection models have BEB"
        print "    - p0, the probability of belonging to the first class of sites with ω<1"
        print "    - p1, the probability of belonging to the second class of sites with ω=1"
        print "    - p2, the probability of belonging to the third class of sites with ω>1)\n"

        if 'BEB' in tree.get_evol_model(model).sites:
            print_selection('BEB', model, tree, verbose)
            
        elif 'NEB' in tree.get_evol_model(model).sites:
            print_selection('NEB', model, tree, verbose)

        print '\n'


def frame_sites(model, tree, models_types, verbose):

    def build_frame(Analysis, model, tree, verbose): 

        plist = {
            True: ['p0','p1','p2'],
            False: ['p2']
        }[verbose]

        for p in plist:
            if p in tree.get_evol_model(model).sites[Analysis]:
                rows = []
                for s in range(len(tree.get_evol_model(model).sites[Analysis]['aa'])):
                    if tree.get_evol_model(model).sites[Analysis][p][s] > 0.95 :
                        rows.append({"Site": tree.get_evol_model(model).sites[Analysis]['aa'][s], "Position": s+1, "p-value": tree.get_evol_model(model).sites[Analysis][p][s]})

                if rows != []:
                    with open('{}_{}.csv'.format(model, p), 'w') as csvfile:
                            writer = csv.DictWriter(csvfile, fieldnames=rows[0].keys())
                            writer.writeheader()
                            for row in rows:
                                writer.writerow(row)

    if model in models_types["site"] :
        if 'BEB' in tree.get_evol_model(model).sites:
            build_frame('BEB', model, tree, verbose)
                
        elif 'NEB' in tree.get_evol_model(model).sites:
            build_frame('NEB', model, tree, verbose)

def details_on_branches(model, tree, models_types):
    if model in models_types["branch"]:
        # branch models have a branches dictionary were keys corresponds to paml_id of nodes in the tree
        # select one of the marked branches
        frg_node = tree.search_nodes(_nid=2)[0]
        frg_pamlid = frg_node.paml_id
        w_frg = bfree.branches[frg_pamlid]['w']
        # select one of the unmarked branches
        bkg_node = tree.search_nodes(_nid=1)[0]
        bkg_pamlid = bkg_node.paml_id
        w_bkg = bfree.branches[bkg_pamlid]['w']
        print 'foreground branches evolving at omega value of %s significantly diferent from %s.' % (w_frg, w_bkg)

def models_types():

    m_types = {"site": ['M1', 'M10', 'M11', 'M12', 'M13', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M8a', 'M9', 'SLR'],
                    "branch": ['b_free', 'b_neut', 'fb'],
                    "branch-site": ['bsA', 'bsA1', 'bsB', 'bsC', 'bsD'],
                    "branch-ancestor": ['fb_anc'],
                    "null": ['M0']}

    return m_types
