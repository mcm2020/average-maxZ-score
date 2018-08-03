#!/usr/bin/env python2

import argparse
from matrix import Matrix
from alignment import Alignment

class Options:

    def __init__(self):
        self.alignment = self.load_args() 

    def load_args(self):
        parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter
        )
        parser.add_argument("-a", "--alignment", required=True,
                            help="\nEnter alignment file.\n\n")

        parser.add_argument("-s", "--symbol", required=True,
                            choices=["protein","nucleotide"],
                            help="\nEnter type of alignment being inputted."
                            "\n\n")
        
        parser.add_argument("-m", "--matrix", required=True,
                            help="\nSelect similarity matrix.\n\n"
                            "| PRESETS: blosum62, blosum90, blosum100,\n"
                            "pam100, pam250, binary.\n\n"
                            "| NOTE: If a preset is not chosen the program\n"
                            "will assume you are loading a file. Make sure\n"
                            "the file you load follows the standard format\n"
                            "set by PAM and BLOSUM.\n\n"
                            "| NOTE: Binary matrices ignore any\n"
                            "similarities among disparate symbols.\n\n"
                            "| NOTE: Nucelotide alignments that are not\n"
                            "converted must use a binary matrix."
                            "\n\n")
        
        parser.add_argument("-d", "--distribution", required=True,
                            help="\nSelect sequences to define the\n" 
                            "background distribution.\n\n"
                            "| PRESETS: self, swiss-prot\n\n"
                            "| NOTE: If a preset is not chosen the program\n"
                            "will assume you are loading a file. Make sure\n"
                            "any file you load follows fasta formatting\n"
                            "and is of the same symbol type (protein or \n"
                            "nucleotide).\n\n"
                            "| NOTE: If self is chosen the sequences\n"
                            "of the inputted alignment file will be used\n"
                            "for the calculation.\n\n")

        parser.add_argument("-c", "--convert", required=True,
                            choices=["yes", "no"],
                            help="\nConvert nucleotide alignment to protein\n" 
                            "alignment?\n\n"
                            "| NOTE: Sequences in alignment must be of equal\n"
                            "length to convert.\n\n"
                            "|NOTE: First open reading frame is used.")
        args = parser.parse_args()

        matrix_input = args.matrix
        if args.symbol == "nucleotide" and args.convert == "no":
            if args.distribution == "swiss-prot":
                exit("Nucleotide alignment cannot use amino acid background "
                     "distribution")
            if matrix_input != "binary":
                exit("Nucleotide alignments that are not converted must use a "
                    "binary matrix")
            matrix_input = "binary-nucleotide"
        else:
            if args.symbol == "protein" and args.convert == "yes":
                exit("Cannot convert protein alignment to nucleotide alignment")
            if matrix_input == "binary":
                matrix_input = "binary-protein"
        sim_matrix = Matrix(matrix_input)

        return Alignment(args.alignment, args.symbol, args.distribution,
                         sim_matrix, args.convert)
    
