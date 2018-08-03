#!/usr/bin/env python2

import numpy as np
from collections import Counter

class ImportanceSampling:

    def __init__(self, alignment, e=.7):
        self.seq_len = alignment.seq_len
        self.sym_set = alignment.sym_set
        self.sym_set_len = alignment.sym_set_len
        self.sim_matrix = alignment.sim_matrix
        self.bg_dist = alignment.bg_dist
        self.num_seqs = alignment.num_seqs
        self.dist_matrix = self.create_dist_matrix(e)
        
    def generate_g0_sample(self):
        # Normalizing bg_dist as probablities commonly do not sum to 1.
        # The sum is checked when bg_dist is created in alignment.py with
        # np.isclose(sum, 1.0)
        # Despite np.isclose() returning true for bg_dist, using random.choice()
        # may still result in a "ValueError: probablities do not sum to 1"
        # This implies that np.isclose() has a higher tolerance than choice().
        # To account for this the vector is normalized.
        bg_dist = np.divide(self.bg_dist, np.sum(self.bg_dist))
        return np.random.choice(self.sym_set, self.num_seqs, p=bg_dist)

    def create_dist_matrix(self, e):
        matrix = np.empty(shape=(self.sym_set_len + 1, self.sym_set_len))
        matrix[0] = np.multiply(np.subtract(1, e),
                                np.divide(self.bg_dist, self.sym_set_len))
        for i in range(self.sym_set_len):
            for j in range(self.sym_set_len):
                var1 = np.divide(self.bg_dist[i], self.sym_set_len)
                var2 = np.multiply(np.subtract(1,e), var1)
                if i==j:
                    matrix[i + 1, j] = np.add(e, var2)
                else:
                    matrix[i + 1, j] = var2
        return matrix
        
    def sym_freqs(self, syms):
        syms_count = Counter(syms)
        freqs = np.empty(shape=(self.sym_set_len))
        for i, sym in enumerate(self.sym_set):
            freqs[i] = syms_count[sym]
        return freqs
    
    def mixture_components(self):
        components = np.empty(shape=(self.sym_set_len + 1))
        components[0] = self.generate_g0_sample
        for i in range(self.sym_set_len):
            components[i+1] = self.generate_g0_sample

    
""" 
    #self.weighted_sym_arr = self.create_weighted_sym_arr()

    # ex: [(0., 0.12), (0.12, .40), (0.40, 73), (0.73, 1.)]
    #          A             C           G           T
    def create_weighted_sym_arr(self):
        ret_arr = np.empty(shape=(self.sym_set_len, 2))
        
        range_end = 0.0
        for i in range(self.sym_set_len):
            range_start = range_end
            range_end = np.add(self.bg_dist[i], range_end)
            ret_arr[i] = (range_start, range_end)
        return ret_arr

    # ran_val should be seeded and within [0, 1)
    def get_aa_weighted(self, ran_val):
        for i, cur_range in enumerate(self.weighted_sym_arr):
            if (np.greater_equal(ran_val, cur_range[0]) and
                np.less(ran_val, cur_range[1])):
                return self.sym_set[i]
        exit("error")
"""
