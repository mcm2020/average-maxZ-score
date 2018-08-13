#!/usr/bin/env python2

import numpy as np
import re
from string import maketrans

# GLOBAL CONSTANTS -------------------------------------------------------------
nucleotide_bases = ['A','C','G','T']
nb_indexed = {nb: i for i, nb in enumerate(nucleotide_bases)}

amino_acids = ['A','R','N','D','C','Q','E','G','H','I',
               'L','K','M','F','P','S','T','W','Y','V',]
aa_indexed = {aa: i for i, aa in enumerate(amino_acids)}

amino_acids_with_ambig = ['A','R','N','D','C','Q','E','G','H','I','L','K',
                          'M','F','P','S','T','W','Y','V','B','Z','X','*',]

# Frequencies generated from the Swiss-Prot protein sequence database as of
# 7-12-18. The db is "high quality, manually annotated and non-redundant."
swiss_prot_bg_dist = np.array(
    [0.08259439291945994, 0.05536917228819641,
     0.04060174708503852, 0.05462634423963772,
     0.013806080104742673, 0.03932612819708233,
     0.06737357357036831, 0.07079896449402255,
     0.022752062190322164, 0.05928018275121051,
     0.09657113918749555, 0.05819965516580402,
     0.024160310383956194, 0.038653182630623874,
     0.04729861257813521, 0.06621618774705094,
     0.05352881699716259, 0.010970035610933597,
     0.029201213397744546, 0.06867040954942194])

# PRETTY PRINTING A FREQUENCY VECTOR IN SORTED (BY FREQUENCY) ORDER ------------
# Useful for debugging.
def print_freq_vector(freqs, symset):
    syms_to_freq = {}
    for i, freq in enumerate(freqs):
        syms_to_freq[symset[i]] = freq
    for k,v in sorted(syms_to_freq.iteritems(), key=lambda x: -x[1]):
        print k, v #"{0:.2f}".format(v*100)
    return ''

# CONVERTING A SEQUENCE OF NUCLEOTIDE BASES TO A PROTEIN USING THE FIRST READING
# FRAME ------------------------------------------------------------------------
pre_tr = 'ATGC'
post_tr = 'TACG'
complement_table = maketrans(pre_tr, post_tr)

codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T',
    'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S',
    'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P',
    'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V',
    'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D',
    'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S',
    'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*',
    'TGG':'W', 'CT':'L', 'GT':'V', 'TC':'S', 'CC':'P', 'AC':'T', 'GC':'A',
    'CG':'R', 'GG':'G', '---':'-', '--':'-', '-':'-',
}

def seq_to_protein(seq):
    seq = seq.rstrip()
    codons = [seq[i:i+3] for i in range(0, len(seq)+1, 3)]
    if codons[-1] == '':
        codons.pop()

    protein = ''
    for codon in codons:
        # 'X' in cases where the last codon is only 1-2 bases and cannot be
        # identified
        protein += codon_table.get(codon, 'X')
    return protein

# READING IN A FASTA FILE AND RETURNING A DICTIONARY OF SEQS MAPPED TO HEADER
# NAMES ------------------------------------------------------------------------

def load_input(filename):
    cur_header = ''
    seqs = {}

    with open(filename) as input:
        for line in input:
            if re.search("^>", line):
                cur_header = re.search(r"^>(\S+)", line).group(1)
                seqs[cur_header] = ''
            else:
                seqs[cur_header] += line.rstrip()

    return seqs

