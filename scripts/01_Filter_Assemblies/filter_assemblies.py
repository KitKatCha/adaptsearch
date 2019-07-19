#!/usr/bin/env python
#coding: utf-8

"""
Usage : python filter_assemblies.py \
            <input,files,comma,separated> \
            <Minimum sequence length> \
            <percent identity> \
            <overlap length> \
            <assembleur>

args[2,3,4] : integers
args[5] : ('trinity', 'velvet_oases')

Author : Eric Fontanillas
last version : 07.2019 by Victor Mataigne
"""

import string, os, sys, re, itertools

"""
This step might actually be skipped, but in case we have to implement the handling
of variations in velvet_oases headers, it will be nice to have the same scheme than
with Trinity

Args:
    - path_in (str) : path/name of the input file (from fasta formatter)
return:
    - - formatted_headers (dict) : {header:sequence} with formatted velvet headers
"""
def format_headers_velvet(path_in):

    formatted_headers = {}

    with open(path_in, 'r') as f_in:
        for header, sequence in itertools.izip_longest(*[f_in]*2):
            h = header.split('_')
            h = '{}_{}_{}'.format(h[1], h[3].split('/')[0], h[5]) # gene_transcript_confidence
            formatted_headers[h] = sequence.rstrip()

    return formatted_headers

"""
Remove redondant transcripts (i.e. transcript from the same locus) from velvet oases ;
criteria 1 - Keep in priority seq with the best "confidence" (parameter in the header)
criteria 2 - If same coverage : choose the longuest sequence (once any "N" have been removed)
Limit of this approach: the transcripts may come from a same locus but may be
not redundant (non-overlapping)

Args:
    - dict_in (dict) : {header:sequence}, content of a Trinity-formatted fasta file
Return :
    - dict_unredundant (dict) : {header:sequence}, i;e input dict without redondant transcripts
"""
def filter_redundancy_velvet(dict_in):

    dict_step_one = {}
    dict_unredundant = {}

    for header, sequence in dict_in.items():
        gene, _, confidence = header.split('_')
        #short_fasta_name = h_split[0]
        countN = sequence.count('N')
        length = len(sequence)
        effective_length = length - countN

        if gene not in list(dict_step_one.keys()):
            dict_step_one[gene] = [[header, sequence, confidence, effective_length]]
        else:
            dict_step_one[gene].append([header, sequence, confidence, effective_length])

    for key in dict_step_one.keys():
        # If one transcript per locus, record directly
        if len(dict_step_one[key]) == 1:
            entry = dict_step_one[key][0]
            name, seq, _, _ = entry
            dict_unredundant[name] = seq

        # If more than one transcript per locus, apply criterias 1 then 2
        elif len(dict_step_one[key]) > 1:
            max_confidence = {} # dict for criteria 1
            max_length = {} # dict for criteria 2

            for entry in dict_step_one[key]:
                name, seq, confidence, effective_length = entry

                max_length[effective_length] = entry
                confidence = float(confidence)

                if confidence not in list(max_confidence.keys()):
                    max_confidence[confidence] = entry
                else:
                    # if several sequences with the same confidence interval : record only the longest one
                    current_seq_length = effective_length
                    yet_recorded_seq_length = max_confidence[confidence][3]
                    if current_seq_length > yet_recorded_seq_length:
                        # Replace the previous recorded entry with the same confidence interval but lower length
                        max_confidence[confidence] = entry

            # Sort keys() for max_confidence dict
            kc = list(max_confidence.keys())
            kc.sort()

            ##Select the best entry
            max_confidence_key = kc[-1]  # criteria 1
            best_entry = max_confidence[max_confidence_key]

            best_fasta_name = best_entry[0]
            best_seq = best_entry[1]
            dict_unredundant[best_fasta_name] = best_seq

    return dict_unredundant

"""
Args:
    - path_in (str) : path/name of the input file
    - trinity (str) : indicates which version of trinity has been used (for using good regex)
Return :
    - formatted_headers (dict) : {header:sequence} with formatted trinity headers
"""
def format_headers_trinity(path_in, trinity=['old', 'last']):

    formatted_headers = {}

    with open(path_in, 'r') as f_in:
        for header, sequence in itertools.izip_longest(*[f_in]*2):
            if trinity == 'old':
                # Catch the cX_gX_iX part of trinity header
                sub_header = re.search('([0-9]+_g[0-9]+_i[0-9]+)', header.rstrip())
                formatted_headers[sub_header.group()] = sequence.rstrip()
            elif trinity == 'last':
                # TODO Catch the DNX_cX_gX_iX part of trinity header
                sub_header = re.search('([0-9]+_c[0-9]+_g[0-9]+_i[0-9]+)', header.rstrip())
                formatted_headers[sub_header.group()] = sequence.rstrip()

    return formatted_headers

"""
Remove redondant transcripts (i.e. transcript from the same locus) from Trinity ;
keeps the longest sequence once any "N" have been removed

Args:
    - dict_in (dict) : {header:sequence}, content of a Trinity-formatted fasta file
Return :
    - dict_unredundant (dict) : {header:sequence}, i;e input dict without redondant transcripts
"""
def filter_redundancy_trinity(dict_in):

    dict_step_one = {}
    dict_unredundant = {}

    for header, sequence in dict_in.items():
        h_split = header.split('_')
        short_fasta_name = h_split[0] # + '_' + h_split[1]
        countN = sequence.count('N')
        length = len(sequence)
        effective_length = length - countN

        if short_fasta_name not in list(dict_step_one.keys()):
            dict_step_one[short_fasta_name] = [[header, sequence, effective_length]]
        else:
            dict_step_one[short_fasta_name].append([header, sequence, effective_length])

    for key in list(dict_step_one.keys()):
        # If one transcript per locus : record directly
        if len(dict_step_one[key]) == 1:
            entry = dict_step_one[key][0]
            name = entry[0]
            seq = entry[1]
            dict_unredundant[name] = seq

        # If more than one transcript per locus : choose the longest sequence (effective_length)
        elif len(dict_step_one[key]) > 1:
            max_length = {}
            for entry in dict_step_one[key]:    ## key = short fasta name    || VALUE = list of list, e.g. :  [[fasta_name1, fasta_seq1],[fasta_name2, fasta_seq2][fasta_name3, fasta_seq3]]
                name = entry[0]
                seq = entry[1]
                effective_length = entry[2]

                ## Bash for [CRITERIA 1]
                max_length[effective_length] = entry

            ## Sort keys() for max_length dict_step_one
            KC = list(max_length.keys())
            KC.sort()

            ## Select the best entry
            max_length_key = KC[-1]  ## [CRITERIA 1]
            best_entry = max_length[max_length_key]

            best_fasta_name = best_entry[0]
            best_seq = best_entry[1]
            dict_unredundant[best_fasta_name] = best_seq

    return dict_unredundant

"""
Args :
    - entry (str) : a RNA sequence
Return :
    - orf (dict) : sequence's open reading frames
"""
def find_orf(entry):
    orf = {}
    orf_length = {}
    stop = ['TAA','TAG','TGA']

    for i in range(0,3):
        pos = i
        orf[i] = [0]

        while pos < len(entry):
            if entry[pos:pos + 3] in stop:
                orf[i].append(pos - 1)
                orf[i].append(pos + 3)
            pos += 3

        orf[i].append(len(entry) - 1)
        orf_length[i] = []

        for u in range(1,len(orf[i])):
            orf_length[i].append(orf[i][u] - orf[i][u - 1] + 1)

        orf[i] = [orf[i][orf_length[i].index(max(orf_length[i]))],orf[i][orf_length[i].index(max(orf_length[i])) + 1]]

    orf_max = {0:max(orf_length[0]),1:max(orf_length[1]),2:max(orf_length[2])}
    orf = orf[max(list(orf_max.keys()), key=(lambda k: orf_max[k]))]

    if orf[0] == 0:
        orf[0] = orf[0] + max(list(orf_max.keys()), key=(lambda k: orf_max[k]))

    return orf

"""
Args :
    - entry (str) : A RNA sequence
Return :
    - seq (str) : reverse complement of the sequence
"""
def reverse_seq(entry):
    nt = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
    seqlist = []

    for i in range(len(entry) -1, -1, -1):
        seqlist.append(nt[entry[i]])

    seq = ''.join(seqlist)

    return seq

"""
Wrapper function for orf finder. Write results in a file instead of returning a
dict because cap3 is called afterwards

Args :
    - dict_in (dict) :
    - file_out (str) : output file path/name
"""
def find_orf_process(dict_in, file_out):

    threshold = 0 #minimal length of the ORF
    f_out = open(file_out, 'w')

    for gene in sorted(list(dict_in.keys())):
        # find longest orf in both strands
        sequence = dict_in[gene]
        high_plus = find_orf(sequence)
        reverse = reverse_seq(sequence)
        high_minus = find_orf(reverse)

        # check threshold for both orf strands
        if high_plus[1] - high_plus[0] > threshold or high_minus[1] - high_minus[0] > threshold:
            # keep the longest one
            if high_plus[1] - high_plus[0] > high_minus[1] - high_minus[0]:
                f_out.write('>{}_{}\n'.format(gene, str(high_plus[1]-high_plus[0]+1)))
                f_out.write(sequence[high_plus[0]:high_plus[1]+1]+'\n')
            else:
                f_out.write('>{}_{}\n'.format(gene, str(high_minus[1]-high_minus[0]+1)))
                f_out.write(sequence[high_minus[0]:high_minus[1]+1]+'\n')

"""
filters the sequences depending on their length after cap3 & makes the sequences
names compatible with the phylogeny workflow

Args :
    - file_in (str) : input file path/name
    - identifier (str) : the species identifier
    - threshold (int) : minimum number of nucleotides for one sequence
    - file_out (str) : output file name
"""
def filter(file_in, identifier, threshold, file_out):

    f_out = open(file_out, "w")

    with open(file_in, "r") as f_in:
        for header, sequence in itertools.izip_longest(*[f_in]*2):
            if len(sequence) - 1 >= threshold :
                name='>{}_{}\n'.format(identifier, header.rstrip().replace('>',''))
                f_out.write("{}".format(name))
                f_out.write("{}".format(sequence)) # no need for '\n' (did not rstrip sequence)

    f_out.close()

def main():
    os.mkdir("outputs")
    script_path = os.path.dirname(sys.argv[0])
    length_seq_max = int(sys.argv[2])
    percent_identity = int(sys.argv[3])
    overlap_length = int(sys.argv[4])

    tmp_names = ['01_', '02_', '03_']

    for name in str.split(sys.argv[1], ','):

        identifier = name.split('.')[0]
        prefix_names = [tmp+identifier for tmp in tmp_names]

        # sequences on single lines
        os.system('cat "{}" | fasta_formatter -w 0 -o "{}"'.format(name, prefix_names[0]))

        # detect the assembly
        assembleur = os.popen('head {} -n 1'.format(prefix_names[0])).read()

        # formatting and filtering isoforms
        if re.search('[0-9]+_g[0-9]+_i[0-9]+', assembleur) != None:
            # Trinity old headers all have the form : 'c[0-9]+_g[0-9]+_i[0-9]+'
            # For the most recent version (07/2019), a second regex must be applied
            if re.search('[0-9]+_c[0-9]+_g[0-9]+_i[0-9]+', assembleur) != None :
                dict_transcripts = format_headers_trinity(prefix_names[0], 'last')
            else:
                dict_transcripts = format_headers_trinity(prefix_names[0], 'old')
            dict_transcripts = filter_redundancy_trinity(dict_transcripts)
        elif re.search('>Locus_[0-9]+_Transcript_[0-9]+/[0-9]+_Confidence_', assembleur) != None:
            dict_transcripts = format_headers_velvet(prefix_names[0])
            dict_transcripts = filter_redundancy_velvet(dict_transcripts)
        else :
            raise ValueError('Wrong assembleur (no Trinity or Velvet Oases)')

        # Pierre guillaume code for keeping the longuest ORF
        find_orf_process(dict_transcripts, prefix_names[1])
        # Apply cap3
        os.system('cap3 {} -p {} -o {}'.format(prefix_names[1], percent_identity, overlap_length))
        # Il faudrait faire un merge des singlets et contigs! TODO
        os.system('zcat -f < "{}.cap.singlets" | fasta_formatter -w 0 -o "{}"'.format(prefix_names[1], prefix_names[2]))
        # Apply pgbrun script filter script TODO length parameter
        filter(prefix_names[2], identifier, length_seq_max, 'outputs/'+name)

    os.mkdir('tmp')
    os.system('mv 01* 02* 03* tmp/')

if __name__ == "__main__":
    main()
