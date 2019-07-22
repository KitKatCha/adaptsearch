#!/usr/bin/env python
# coding: utf-8
# Author : Victor Mataigne

import argparse, itertools, pickle

def list_with_max_score(list_of_hits):
    """ Among a list of blast hits of the same query, returns the one which has the highest score """
    max_score = 0
    ind_max_score = 0
    i = 0

    for hit in list_of_hits:
        if float(hit[11]) > max_score:
            max_score, ind_max_score = hit[11], i
        i += 1

    return list_of_hits[ind_max_score]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('matches', help='diamond blastp output file (tabular) between two species')
    parser.add_argument('nucleic_file_db', help='Sequence file used as DB for the first blast')
    parser.add_argument('nucleic_file_q', help='Sequence used as query for the first blast')
    parser.add_argument('file_subname', help='keyword for output file name')
    parser.add_argument('step', help='1 or 2 according to which blast has been performed')
    parser.add_argument('method', choices=['tblastx', 'diamond'], help='alignment tool (tblastx or diamond)')
    args = parser.parse_args()

    print 'Keeping best hits in {}'.format(args.matches)

    # read tab file in a list (a line = a list elem)
    list_hits = []
    dic_hits_common = {} # unique matches for query
    dic_hits_common_db = {} # unique matches for db
    with open(args.matches, 'r') as hits:
        for hit in hits.readlines():
            h = hit.strip('\n')
            h2 = h.split('\t')
            list_hits.append(h2)

            if args.method == 'diamond':
                if h2[0][0:-6] not in dic_hits_common.keys():
                        dic_hits_common[h2[0][0:-6]] = []
                if h2[1][0:-6] not in dic_hits_common_db.keys():
                        dic_hits_common_db[h2[1][0:-6]] = []

            elif args.method == 'tblastx':
                if h2[0] not in dic_hits_common.keys():
                    dic_hits_common[h2[0]] = []
                if h2[1] not in dic_hits_common_db.keys():
                    dic_hits_common_db[h2[1]] = []

    # Gather in a list of lists elems with common query
    for hit in list_hits:
        if args.method == 'diamond':
            dic_hits_common[hit[0][0:-6]].append(hit)
            dic_hits_common_db[hit[1][0:-6]].append(hit)

        elif args.method == 'tblastx':
            dic_hits_common[hit[0]].append(hit)
            dic_hits_common_db[hit[1]].append(hit)

    # Keep only the best hits in queries
    list_best_hits_q = []
    for list_hits in dic_hits_common.values():
        list_best_hits_q.append(list_with_max_score(list_hits))

    # Keep only the best hits in db
    list_best_hits_db = []
    for list_hits in dic_hits_common_db.values():
        list_best_hits_db.append(list_with_max_score(list_hits))

    del list_hits

    # This dict (exported then with pickle) stores the best hit in the db for the query
    # A similar dict is built after the second blast, in which query and db are switched
    # The comparison of the dicts allow to spot RBH
    dico_best_hits_q = {}
    for hit in list_best_hits_q:
        if args.method == 'diamond':
            dico_best_hits_q[hit[0][0:-6]] = hit[1][0:-6]
        elif args.method == 'tblastx':
            dico_best_hits_q[hit[0]] = hit[1]

    n = 'dict_best_hits_from_blast_{}'.format(args.step)
    pickle_dic_besthits_q = open(n, 'w')
    pickle.dump(dico_best_hits_q, pickle_dic_besthits_q)
    pickle_dic_besthits_q.close()

    ## Other temp files :

    # Make big dictionary with initial query fasta file
    initial_seqs_q = {}
    with open(args.nucleic_file_q, 'r') as nf :
        for h, s in itertools.izip_longest(*[nf]*2):
            header = h.strip('\n')
            sequence = s.strip('\n')
            initial_seqs_q[header] = sequence

    # Make big dictionary with initial DB fasta file
    initial_seqs_db = {}
    with open(args.nucleic_file_db, 'r') as nf :
        for h, s in itertools.izip_longest(*[nf]*2):
            header = h.strip('\n')
            sequence = s.strip('\n')
            initial_seqs_db[header] = sequence

    # Write best_hits from query with nucleic sequence in output file
    name = 'best_hits_q_blast{}_{}'.format(args.step, args.file_subname)
    output = open(name, 'w')
    for hit in list_best_hits_q:
        if args.method == 'diamond':
            output.write('>'+hit[0][0:-6]+'\n')
            output.write(initial_seqs_q['>'+hit[0][0:-6]]+'\n')
        elif args.method == 'tblastx':
            output.write('>'+hit[0]+'\n')
            output.write(initial_seqs_q['>'+hit[0]]+'\n')
    output.close()

    # Write best_hits on db with nucleic sequence in output file
    name = 'best_hits_db_blast{}_{}'.format(args.step, args.file_subname)
    output = open(name, 'w')
    for hit in list_best_hits_db:
        if args.method == 'diamond':
            output.write('>'+hit[1][0:-6]+'\n')
            output.write(initial_seqs_db['>'+hit[1][0:-6]]+'\n')
        elif args.method == 'tblastx':
            output.write('>'+hit[1]+'\n')
            output.write(initial_seqs_db['>'+hit[1]]+'\n')

    output.close()

    del initial_seqs_q
    del initial_seqs_db


if __name__ == "__main__":
    main()