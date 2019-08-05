#!/usr/bin/env python
# coding: utf-8
# Author : Victor Mataigne

import argparse, pickle, itertools

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('besthits_file1', help='')
    parser.add_argument('besthits_file2', help='')
    parser.add_argument('out_name_id', help='')
    args = parser.parse_args()

    # Open dict of best hits
    file_best_hit_dict_q = open('dict_best_hits_from_blast_1')
    file_best_hit_dict_db = open('dict_best_hits_from_blast_2')
    best_hit_dict_q = pickle.load(file_best_hit_dict_q)
    best_hit_dict_db = pickle.load(file_best_hit_dict_db)
    file_best_hit_dict_q.close()
    file_best_hit_dict_db.close()

    best_h1 = {}
    with open(args.besthits_file1, 'r') as bh1 :
        for h, s in itertools.izip_longest(*[bh1]*2):
            header = h.strip('>\n')
            sequence = s.strip('\n')
            best_h1[header] = sequence

    best_h2 = {}
    with open(args.besthits_file2, 'r') as bh2 :
        for h, s in itertools.izip_longest(*[bh2]*2):
            header = h.strip('>\n')
            sequence = s.strip('\n')
            best_h2[header] = sequence

    # Find RBH:
    reverse_best_hit_dict_db = dict((v,k) for k,v in best_hit_dict_db.iteritems())

    rbh = set(best_hit_dict_q.items()).intersection(set(reverse_best_hit_dict_db.items()))

    if len(rbh) > 0:

        s = args.besthits_file1.split('_')
        suffix = s[4] + '_' + s[5]
        out_name = 'RBH_{}_dna.fasta'.format(args.out_name_id)
        output = open(out_name, 'w')

        for pairwise_couple in rbh :
            output.write('>'+pairwise_couple[0]+'\n')
            output.write(best_h1[pairwise_couple[0]]+'\n')
            output.write('>'+pairwise_couple[1]+'\n')
            output.write(best_h2[pairwise_couple[1]]+'\n')

        output.close()

if __name__ == "__main__":
    main()