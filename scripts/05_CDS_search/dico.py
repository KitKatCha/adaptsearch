import string, itertools

def dico(file_in):
    dicoco = {}

    with open(file_in, "r") as f_in:
        for name, query in itertools.izip_longest(*[f_in]*2):
            fasta_name_query = name.rstrip()
            fasta_seq_query = query.rstrip()
            dicoco[fasta_name_query] = fasta_seq_query

    return(dicoco)
