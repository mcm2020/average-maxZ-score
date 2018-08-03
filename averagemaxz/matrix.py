#!/usr/bin/env python2

import config
import matrices
import numpy as np

class Matrix:

    def __init__(self, filename):
        if filename in stored_matrices:
            self.matrix = stored_matrices[filename]
        else:
            self.matrix = self.load_matrix(filename)
    
    def load_matrix(self, filename):
        line_entries = []
        with open(filename) as input:
            for line in input:
                if line.startswith(("#", " ", "B", "Z", "X", "*")) == False :
                    line = line.split()

                    # column 1 is amino acid ID
                    # last 4 columns refer to ambiguous amino acids
                    line = line[1:-4] 

                    line = [int(x) for x in line]
                    line_entries.append(line)

        return np.array(line_entries)
        
stored_matrices = {"binary-protein": matrices.binary_protein,
                   "binary-nucleotide": matrices.binary_nucleotide,
                   "blosum62": matrices.blosum62,
                   "blosum90": matrices.blosum90,
                   "blosum100": matrices.blosum100,
                   "pam100": matrices.pam100,
                   "pam250": matrices.pam250}


