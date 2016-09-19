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
    `CODON_TO_AA` (`numpy.ndarray` of int, lenght `N_CODON`)
        Element `x` gives index for amino acid encoded by `x`.
    `PURINES` (`frozenset`)
        Set of one-letter nucleotides for purines
    `PYRIMIDINES` (`frozenset`)
        Set of one-letter nucleotides for pyrimidines
    `CODON_TRANSITION` (`numpy.ndarray` of bools, shape `(N_CODON, NCODON)`)
        Element `[x][y]` is `True` iff `x` and `y` differ by 1 transition.
    `CODON_SINGLEMUT` (`numpy.ndarray` of bools, shape (`N_CODON, N_CODON)`)
        Element `[x][y]` is `True` iff `x` and `y` differ by 1 nt mutation.
    `CODON_NT_MUT` (`numpy.ndarray` of bools, shape `(N_NT, N_CODON, N_CODON)`)
        Element `[w][x][y]` is `True` iff `x` converts to `y` by single 
        nucleotide mutation `w`.
    `CODON_NT` (`numpy.ndarray` of bools, shape `(3, N_NT, N_CODON)`)
        Element `[j][w][x]` is `True` iff nt `j` of codon `x` is `w`.
    `CODON_NONSYN` (`numpy.ndarray` of bools, shape `(N_CODON, N_CODON)`)
        Element `[x][y]` is `True` iff mutating  `x` to `y` is nonsynonymous
"""


import re
import inspect
import scipy
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
CODON_TO_AA = []
i = 0
for nt1 in NT_TO_INDEX.keys():
    for nt2 in NT_TO_INDEX.keys():
        for nt3 in NT_TO_INDEX.keys():
            codon = nt1 + nt2 + nt3
            aa = str(Bio.Seq.Seq(codon).translate())
            if aa != '*':
                CODON_TO_INDEX[codon] = i
                INDEX_TO_CODON[i] = codon
                CODON_TO_AA.append(AA_TO_INDEX[aa])
                i += 1
N_CODON = len(CODON_TO_INDEX)
CODON_TO_AA = scipy.array(CODON_TO_AA, dtype='int')
assert len(CODON_TO_INDEX) == len(INDEX_TO_CODON) == len(CODON_TO_AA) == N_CODON

PURINES = frozenset(['A', 'G'])
PYRIMIDINES = frozenset(['C', 'T'])
assert PURINES.union(PYRIMIDINES) == frozenset(NT_TO_INDEX.keys())

CODON_TRANSITION = scipy.full((N_CODON, N_CODON), False, dtype='bool')
CODON_SINGLEMUT = scipy.full((N_CODON, N_CODON), False, dtype='bool')
CODON_NT_MUT = scipy.full((N_NT, N_CODON, N_CODON), False, dtype='bool')
CODON_NT = scipy.full((3, N_NT, N_CODON), False, dtype='bool')
CODON_NONSYN = scipy.full((N_CODON, N_CODON), False, dtype='bool')
for (x, codonx) in INDEX_TO_CODON.items():
    for (i, ntx) in enumerate(codonx):
        for w in range(N_NT):
            if INDEX_TO_NT[w] == ntx:
                CODON_NT[i][w][x] = True
    for (y, codony) in INDEX_TO_CODON.items():
        if CODON_TO_AA[x] != CODON_TO_AA[y]:
            CODON_NONSYN[x][y] = True
        diffs = [(ntx, nty) for (ntx, nty) in zip(codonx, codony) if ntx != nty]
        if len(diffs) == 1:
            (ntx, nty) = diffs[0]
            CODON_SINGLEMUT[x][y] = True
            CODON_NT_MUT[NT_TO_INDEX[nty]][x][y] = True
            if ((ntx in PURINES and nty in PURINES) or 
                    ((ntx in PYRIMIDINES and nty in PYRIMIDINES))):
                CODON_TRANSITION[x][y] = True
assert CODON_SINGLEMUT.sum() == CODON_NT_MUT.sum() > CODON_TRANSITION.sum()
assert CODON_NT.sum() == N_CODON * 3

# delete variables so they aren't in namespace if import * used on this module
del i, nt1, nt2, nt3, codon, codonx, codony, x, y, ntx, nty, w, diffs
# following lines needed because list comprehension indices remain
# in Python2 but not Python3, and we want to be compatible with both
if 'nt' in locals():
    del nt
if 'aa' in locals():
    del aa

# make sure we deleted all variables not all upper case
for (key, value) in list(locals().items()):
    if not (inspect.ismodule(value) or re.search('__.+__', key)):
        # is a user-defined variable if we made it here
        assert key.isupper(), "Failed to delete variable: {0}".format(key)
