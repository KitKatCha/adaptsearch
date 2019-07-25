#!/usr/bin/env python
# coding: utf-8
# Author : Victor Mataigne

import itertools, argparse, os

"""
Usage : ./S01_run_first_blast.py <input,files,comma,separated> <evalue> <[tblastx, diamond]>]

WARNING : input files must be at the path than the script

IMPROVMENTS :
    - Maybe a bit of code factoring
    - See if it possible to avoid build several times the same db
"""

# The script (and S03_run_second_blast.py as well) must be launched with the python '-W ignore' option if tested with planemo

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('files', help='fasta files separated by commas (MUST be in script directory')
    parser.add_argument('evalue', help='evalue for blast')
    parser.add_argument('method', choices=['tblastx', 'diamond'], help='alignment tool (tblastx or diamond)')
    args = parser.parse_args()

    in_files = args.files.split(',')

    if args.method == 'diamond':
        os.mkdir('translated_seqs')
        in_files_translated = []
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC

        # From every sequence, make three sequences (translations in the three reading frames)
        print 'Translating every sequence in all reading frames ...'
        for file in in_files:
            #basename_f = file.split('/')[-1]
            name = 'translated_seqs/%s' %file
            in_files_translated.append(name)
            translated_file = open(name, 'w')
            with open(file, 'r') as f:
                for gene, seq in itertools.izip_longest(*[f]*2):
                    s = Seq(seq.strip('\n').upper(), IUPAC.ambiguous_dna)
                    translated_file.write(gene.strip('\n')+'_orf_1\n')
                    translated_file.write(s.translate()._data+'\n')
                    translated_file.write(gene.strip('\n')+'_orf_2\n')
                    translated_file.write(s[1:].translate()._data+'\n')
                    translated_file.write(gene.strip('\n')+'_orf_3\n')
                    translated_file.write(s[2:].translate()._data+'\n')
            translated_file.close()

        # Make the list of all pairwise combinations
        list_pairwise = itertools.combinations(in_files_translated, 2)

    elif args.method == 'tblastx':
        list_pairwise = itertools.combinations(in_files, 2)

    else :
        print 'Mispecified alignment tool'
        exit()

    os.mkdir('outputs_RBH_dna')

    # Main loop

    if args.method == 'diamond':
        for pairwise in list_pairwise:
            print "\n================\nPair of species:"
            print pairwise

            sp1, sp2 = pairwise[0].split('/')[1].split('.')[0], pairwise[1].split('/')[1].split('.')[0]
            sub_directory_name = sp1 + '_' + sp2

            os.mkdir('blast_%s' %sub_directory_name)

            print 'Running first blast with Diamond ...'

            # # Run diamond
            os.system('diamond makedb --in %s -d %s --quiet >> log_diamond.log' %(pairwise[1], sp2))
            os.system('diamond blastp -q %s -d %s --max-target-seqs 1 -o matches_blast1_%s -e %s --more-sensitive --quiet >> log_diamond.log' %(pairwise[0], sp2, sub_directory_name, args.evalue))

            # # tabular output :
            # # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

            a = pairwise[1].replace('translated_seqs/', '')
            b = pairwise[0].replace('translated_seqs/', '')

            # # There is a chance to have no hits returned
            if os.path.getsize('matches_blast1_%s' %sub_directory_name) == 0:
                print 'No hits found. Processing next species pair ...'
            else :

                # Record only one best_hit per transcript (best of the 6 orfs)
                os.system('python S02_04_keep_one_hit_from_blast.py matches_blast1_%s %s %s %s %s %s' %(sub_directory_name, a, b, sub_directory_name, '1', args.method))

                # 2d blast with only best hits as db
                print 'Running second blast with Diamond ... '
                print sub_directory_name, pairwise[0]
                os.system('python -W ignore S03_run_second_blast.py best_hits_db_blast1_%s %s %s %s %s' %(sub_directory_name, pairwise[0], sub_directory_name, args.evalue, args.method))

                # Record only one best_hit per transcript (best of the 6 orfs)
                os.system('python S02_04_keep_one_hit_from_blast.py matches_blast2_%s %s %s %s %s %s' %(sub_directory_name, b, a, sub_directory_name, '2', args.method))

                # Find Reciprocical Best Hits
                name1 = 'best_hits_q_blast1_{}'.format(sub_directory_name)
                name2 = 'best_hits_q_blast2_{}'.format(sub_directory_name)
                os.system('python S05_find_rbh.py %s %s %s' %(name1, name2, sub_directory_name))

            os.system('mv log_diamond.log blast_%s' %sub_directory_name)
            os.system('rm -f *.dmnd')

            # Those files exist obly if hits were found during the first blast
            if os.path.getsize('matches_blast1_%s' %sub_directory_name) != 0:
                os.system('mv *best_hits* ./blast_%s' %sub_directory_name)
                os.system('mv RBH* outputs_RBH_dna')

            os.system('mv matches_blast* ./blast_%s' %(sub_directory_name))

        #os.mkdir('translated_seqs')
        #os.system('mv translated*.fasta ./translated_seqs')

    elif args.method == 'tblastx':
        for pairwise in list_pairwise:
            print "\n================\nPair of species:"
            print pairwise

            sp1, sp2 = pairwise[0].split('.')[0], pairwise[1].split('.')[0]
            sub_directory_name = sp1 + '_' + sp2
            os.mkdir('blast_%s' %sub_directory_name)

            print 'Running first tblastx ...'

            # Run diamond
            os.system('formatdb -i %s -p F -o T >> log_tblastx.log' %(pairwise[1]))
            os.system('blastall -p tblastx -d %s -i %s -o matches_blast1_%s -T F -e %s -F "mS" -b1 -v1 -K 1 -m 8 >> log_tblastx.log' %(pairwise[1], pairwise[0], sub_directory_name, args.evalue))

            # tabular output :
            # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

            # There is a chance to have no hits returned
            if os.path.getsize('matches_blast1_%s' %sub_directory_name) == 0:
                 print 'No hits found. Processing next species pair ...'
            else:

                # Record only one best_hit per transcript (best of the 6 orfs)
                os.system('python S02_04_keep_one_hit_from_blast.py matches_blast1_%s %s %s %s %s %s' %(sub_directory_name, pairwise[1], pairwise[0], sub_directory_name, '1', args.method))

                # 2d blast with only best hits as db
                print 'Running second blast with Diamond ... '

                os.system('python S03_run_second_blast.py best_hits_db_blast1_%s %s %s %s %s' %(sub_directory_name, pairwise[0], sub_directory_name, args.evalue, args.method))

                # Record only one best_hit per transcript (best of the 6 orfs)
                os.system('python S02_04_keep_one_hit_from_blast.py matches_blast2_%s %s %s %s %s %s' %(sub_directory_name, pairwise[0], pairwise[1], sub_directory_name, '2', args.method))

                # Find Reciprocical Best Hits
                name1 = 'best_hits_q_blast1_{}'.format(sub_directory_name)
                name2 = 'best_hits_q_blast2_{}'.format(sub_directory_name)
                os.system('python S05_find_rbh.py %s %s %s ' %(name1, name2, sub_directory_name))

            os.system('mv log_tblastx.log ./blast_%s' %sub_directory_name)
            os.system('rm -f *.nhr')
            os.system('rm -f *.nin')
            os.system('rm -f *.nsd')
            os.system('rm -f *.nsi')
            os.system('rm -f *.nsq')

            # Those files exist obly if hits were found during the first blast
            if os.path.getsize('matches_blast1_%s' %sub_directory_name) != 0:
                os.system('mv *best_hits* ./blast_%s' %sub_directory_name)
                os.system('mv RBH* outputs_RBH_dna')

            os.system('mv matches_blast* ./blast_%s' %(sub_directory_name))

if __name__ == "__main__":
    main()
