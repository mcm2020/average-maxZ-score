# Average MaxZ Score Calculator
**Input:** Multiple Sequence Alignment (MSA) in FASTA file format.

**Output:** Average maxZ score of the MSA.

A statistical measure for gauging the relative polymorphism of sequences in a multiple sequence alignment. Useful for the rank ordering of a defined group of MSA's. 

The maximum Z-score is [a statistical measure for positional conservation (Ahola, Aittokallio, Vinhinen, & Uusipaikka, 2006)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-7-484).

*Ahola, V., Aittokallio, T., Vihinen, M., & Uusipaikka, E. (2006). A statistical score for assessing the quality of multiple sequence alignments. BMC Bioinformatics, 7(1), 484. doi:10.1186/1471-2105-7-484*

## Installation
### To Install: 
```pip install .```

### To Uninstall: 
```pip uninstall averagemaxz```

## Usage
```
maxz [-h/--help] -a/--alignment ALIGNMENT -s {protein,nucleotide} -m/--matrix MATRIX -d/--distribution DISTRIBUTION -c/--convert {yes,no}
```

### Arguments
#### -h/--help *(optional)*
- show help message and exit.

#### -a/--alignment ALIGNMENT
- Enter alignment file.

- NOTE: Make sure file follows FASTA formatting.

#### -s/--symbol {protein,nucleotide}
- Enter type of alignment being inputted.

#### -m/--matrix MATRIX
- Select similarity matrix.

- PRESETS: blosum62, blosum90, blosum100, pam100, pam250, binary.

- NOTE: If a preset is not chosen the program will assume you are loading a file. Make sure the file you load follows the standard format set by PAM and BLOSUM.

- NOTE: Binary matrices ignore any similarities among disparate symbols.

- NOTE: Nucelotide alignments that are not converted must use a binary matrix.

#### -d/--distribution DISTRIBUTION
- Select sequences to define the background distribution.

- PRESETS: self, swiss-prot

- NOTE: If a preset is not chosen the program will assume you are loading a file. Make sure any file you load follows FASTA formatting and is of the same symbol type (protein or nucleotide).

- NOTE: If self is chosen the sequences of the inputted alignment file will be used

#### -c/--convert {yes,no}
- Convert nucleotide alignment to protein alignment?

- NOTE: Sequences in alignment must be of equal length to convert.

- NOTE: First open reading frame is used for conversion.
