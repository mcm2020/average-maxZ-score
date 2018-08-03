#!/usr/bin/env python2

import config
from collections import Counter
import numpy as np

class MaxZ:
    
    def __init__(self, alignment, pos):
        self.observed_syms = [seq[pos] for seq in alignment.seqs.itervalues()
                              if seq[pos] != '-']       
        self.num_observed_syms = len(self.observed_syms)
        self.sym_set = alignment.sym_set
        self.sym_set_len = len(self.sym_set)
        self.relative_freqs = self.calc_relative_freqs()
        self.bg_dist = alignment.bg_dist
        self.sim_matrix = alignment.sim_matrix
        self.covariance_matrix = self.calc_covariance_matrix_h0()
        self.calc_maxz() # self.score, self.sym

    def calc_relative_freqs(self):
        relative_freqs = np.empty(shape=(self.sym_set_len))
        syms_count = Counter(self.observed_syms)
        for i, sym in enumerate(self.sym_set):
            relative_freqs[i] = np.true_divide(syms_count[sym],
                                               self.num_observed_syms)
        return relative_freqs

    def calc_covariance_matrix_h0(self):
        b = self.relative_freqs
        beta0 = self.bg_dist
        matrix = np.empty(shape=(self.sym_set_len, self.sym_set_len)) 
        for i in range(self.sym_set_len):
            for j in range(self.sym_set_len):
                if i == j:
                    kronecker = 1
                else:
                    kronecker = 0;
                sym_prob1 = beta0[i]
                sym_prob2 = beta0[j]
                numerator = np.multiply(sym_prob1, np.subtract(kronecker,
                                                               sym_prob2))
                entry = np.divide(numerator, self.num_observed_syms)
                matrix[i, j] = entry
        return matrix

    # helper for calc_maxz()
    def calc_z_score_of_sym(self, i):
        b = self.relative_freqs
        beta0 = self.bg_dist
        sim_matrix_row = self.sim_matrix[i]
        cov_matrix = self.covariance_matrix
        numerator = np.dot(sim_matrix_row, np.subtract(b, beta0))
        variance = np.dot(np.dot(sim_matrix_row, cov_matrix), sim_matrix_row)
        std_deviation = np.sqrt(variance)
        return np.true_divide(numerator, std_deviation)
    
    def calc_maxz(self):
        z_scores = np.empty(shape=(self.sym_set_len))
        for i in range(self.sym_set_len):
            # Can't calculate z_score of symbol if symbol has 0 presence in the
            # background distribution
            if self.bg_dist[i] == 0:
                z_scores[i] = np.NINF # negative infinity
            else:
                z_scores[i] = self.calc_z_score_of_sym(i)
        self.score = np.max(z_scores)
        self.sym = self.sym_set[z_scores.argmax()]
        return z_scores
    
