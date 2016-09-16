"""Defines constants and nucleotide / amino acid / codon indices.

Constants defined:
    `ALMOST_ZERO` (float)
        lower limit for values that should remain > 0
    `N_NT` (int)
        number of nucleotides
    `N_AA` (int)
        number of amino acids
    `N_CODON` (int) 
        number of codons
    `NT_TO_INDEX` (dict)
        mapping of one-letter nucleotides to integer indices
    `INDEX_TO_NT` (dict)
        inverse of `NT_TO_INDEX`
    `AA_TO_INDEX` (dict)
        mapping of one-letter amino acids to integer indices
    `INDEX_TO_AA` (dict)
        inverse of `AA_TO_INDEX`
    `CODON_TO_INDEX` (dict):
        mapping of non-stop codons to integer indices
    `INDEX_TO_CODON` (dict)
        inverse of `CODON_TO_INDEX`
    `TRANSLATE_BY_INDEX` (dict)
        mapping of codon indices to amino-acid indices
"""

import Bio.Alphabet.IUPAC
import Bio.Seq


ALMOST_ZERO = 1.0e-5

INDEX_TO_NT = dict(enumerate(sorted(Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters)))
NT_TO_INDEX = dict([(nt, i) for (i, nt) in INDEX_TO_NT.items()])
N_NT = len(INDEX_TO_NT)
assert len(INDEX_TO_NT) == len(NT_TO_INDEX) == N_NT

INDEX_TO_AA = dict(enumerate(sorted(Bio.Alphabet.IUPAC.IUPACProtein.letters)))
AA_TO_INDEX = dict([(aa, i) for (i, aa) in INDEX_TO_AA.items()])
N_AA = len(INDEX_TO_AA)
assert len(INDEX_TO_AA) == len(AA_TO_INDEX) == N_AA

CODON_TO_INDEX = {}
INDEX_TO_CODON = {}
TRANSLATE_BY_INDEX = {}
i = 0
for nt1 in NT_TO_INDEX.keys():
    for nt2 in NT_TO_INDEX.keys():
        for nt3 in NT_TO_INDEX.keys():
            codon = nt1 + nt2 + nt3
            aa = str(Bio.Seq.Seq(codon).translate())
            if aa != '*':
                CODON_TO_INDEX[codon] = i
                INDEX_TO_CODON[i] = codon
                TRANSLATE_BY_INDEX[i] = AA_TO_INDEX[aa]
                i += 1
N_CODON = len(CODON_TO_INDEX)
assert len(CODON_TO_INDEX) == len(INDEX_TO_CODON) == len(TRANSLATE_BY_INDEX) == N_CODON

# delete variables so they aren't in namespace if import * used
del i, nt1, nt2, nt3, codon
# following lines needed because list comprehension indices remain
# in Python2 but not Python3, and we want to be compatible with both
if 'nt' in locals():
    del nt
if 'aa' in locals():
    del aa
