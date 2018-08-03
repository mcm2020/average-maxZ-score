#!/usr/bin/env python2

import config
from maxz import MaxZ
import numpy as np
from collections import Counter
import re

class Alignment:
    
    def __init__(self, aligned_file, sym_set, bg_dist, sim_matrix, convert):
        self.seqs = config.load_input(aligned_file)

        if sym_set == "nucleotide":
            if convert == "no":
                self.sym_set = config.nucleotide_bases
            else:
                self.nucleotides_to_amino_acids()
                self.sym_set = config.amino_acids
        else:
            self.sym_set = config.amino_acids
            
        self.seq_len = self.calc_seq_len()
        self.sym_set_len = len(self.sym_set)
        self.num_seqs = len(self.seqs)
        self.sim_matrix = sim_matrix.matrix

        self.preset_bg_dists = {"self": self.calc_bg_dist_aggregate_avg(),
                                "swiss-prot": config.swiss_prot_bg_dist}
        
        if bg_dist in self.preset_bg_dists:
            self.bg_dist = self.preset_bg_dists[bg_dist]
        else:
            self.bg_dist = self.calc_bg_dist_aggregate_avg(bg_dist)

    def calc_seq_len(self):
        length = len(self.seqs.values()[0])
        if any(len(seq) != length for seq in self.seqs.itervalues()):
            exit("sequences must be of identical length")
        return length

    def nucleotides_to_amino_acids(self):
        for header, seq in self.seqs.iteritems():
            self.seqs[header] = config.seq_to_protein(seq)
            
        # checking for unequal seq lengths in converted alignment
        max_len = len(max(self.seqs.itervalues(), key=(lambda x: len(x))))
        for seq in self.seqs.itervalues():
            if len(seq) != max_len:
                exit("error converting to protein sequence")           

    # Calculates background distribution by acquiring aggregate symbol
    # frequencies across all sequences, and then dividing by the size of the
    # symbol alphabet. Preferable to the method below, as this method only
    # performs 20 floating point operations.
    def calc_bg_dist_aggregate_avg(self, filename=None):
        if filename is None:
            seqs = self.seqs
        else:
            seqs = config.load_input(filename)
            
        sym_counter = Counter()
        for seq in seqs.itervalues():
            sym_counter.update(seq) # counting symbols in each seq
                
        if self.sym_set == config.amino_acids:
            # these ambigous amino acids will not be included in the dist.
            ambiguous_amino_acids = ['B', 'Z', 'X', '*', '-']            
            num_ambig_amino_acids = 0
            for sym in ambiguous_amino_acids:
                num_ambig_amino_acids += sym_counter[sym]
        else:
            num_ambig_amino_acids = sym_counter['-']
            
        total_num_syms = sum(sym_counter.values()) - num_ambig_amino_acids
            
        freqs = np.empty(shape=(self.sym_set_len))
        freqs_sum = 0.0 #should add up to 1.0
        for i, sym in enumerate(self.sym_set):
            freqs[i] = np.true_divide(sym_counter[sym], total_num_syms)
            freqs_sum = np.add(freqs_sum, freqs[i])

        if not np.isclose(freqs_sum, 1.0):
            exit("background distribution frequencies do not add up to one")
        return freqs

    def avg_maxz_score(self):
        scores = []
        for i in range(self.seq_len):
            # if pos is not all gaps
            if self.check_num_observed_syms(i) > 0:
                cur_maxz = MaxZ(self, i)
                scores.append(cur_maxz.score)
        return np.mean(scores)

    def check_num_observed_syms(self, pos):
        count = 0
        for seq in self.seqs.itervalues():
            if seq[pos] != '-':
                count += 1
        return count
"""
    # Calculates background distribution using the following method:
    # Calculates symbol freqs within each sequence (accounting for gaps), and
    # then averages the symbol freqs of each sequence together.
    # Note, this method is not used as it involves many floating point
    # operations, which is an issue considering floating point imprecision.
    def calc_bg_dist_avg_of_freqs(filename):
        file_seqs = config.load_input(filename)
        freqs_each_seq = np.empty(shape=(self.sym_set_len, len(file_seqs)))
        for j, seq in enumerate(file_seqs.itervalues()):
            seq_counter = Counter(seq)
            num_gaps = seq_counter["-"]
            seq_len_no_gaps = np.subtract(len(seq), num_gaps)
            for i, sym in enumerate(self.sym_set):
                sym_freq = np.true_divide(seq_counter[sym], seq_len_no_gaps)
                freqs_each_seq[i, j] = sym_freq
        freqs = np.empty(shape=(self.sym_set_len))
        total = 0.0
        for i, seqs in enumerate(freqs_each_seq):
            freqs[i] = np.mean(freqs_each_seq[i])
            total += freqs[i]
        if total != 1.0:
            sys.exit("error")
        return freqs
"""
