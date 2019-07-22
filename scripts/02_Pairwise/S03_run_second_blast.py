#!/usr/bin/env python
# coding: utf-8
# Author : Victor Mataigne

import itertools, argparse, os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('query_file', help='fasta file (to translate) for query')
    parser.add_argument('db_file', help='fasta files (already translated) for db')
    parser.add_argument('file_subname', help='keyword for output file name')
    parser.add_argument('evalue', help='evalue for blast')
    parser.add_argument('method', choices=['tblastx', 'diamond'], help='alignment tool (tblastx or diamond)')
    args = parser.parse_args()

    if args.method == 'diamond':
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        # Traduire les best hits
        f_name = 'translated_seqs/{}'.format(args.query_file)
        translated_file = open(f_name, 'w')
        with open(args.query_file, 'r') as file:
            for gene, seq in itertools.izip_longest(*[file]*2):
                s = Seq(seq.strip('\n').upper(), IUPAC.ambiguous_dna)
                translated_file.write(gene.strip('\n')+'_orf_1\n')
                translated_file.write(s.translate()._data+'\n')
                translated_file.write(gene.strip('\n')+'_orf_2\n')
                translated_file.write(s[1:].translate()._data+'\n')
                translated_file.write(gene.strip('\n')+'_orf_3\n')
                translated_file.write(s[2:].translate()._data+'\n')
        translated_file.close()

        os.system('diamond makedb --in %s -d %s --quiet >> log_diamond.log' %(args.db_file, args.db_file.split('/')[1]))
        os.system('diamond blastp -q %s -d %s --max-target-seqs 1 -o matches_blast2_%s -e %s --more-sensitive --quiet >> log_diamond.log' %(f_name, args.db_file.split('/')[1], args.file_subname, args.evalue))

    elif args.method == 'tblastx':
        os.system('formatdb -i %s -p F -o T >> log_tblastx.log' %(args.db_file))
        os.system('blastall -p tblastx -d %s -i %s -o matches_blast2_%s -T F -e %s -F "mS" -b1 -v1 -K 1 -m 8 >> log_tblastx.log' %(args.db_file, args.query_file, args.file_subname, args.evalue))

    else :
        print 'Mispecified alignment tool'
        exit()

if __name__ == "__main__":
    main()
