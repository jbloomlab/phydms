"""Defines constants and nucleotide / amino acid / codon indices.

*ALMOST_ZERO* : used as a lower limit for values that should remain > 0.

*ALMOST_ONE* : used as an upper limit for values that should remain > 1.

*LARGE_NUMBER* : used as an upper limit for values that can be large.

*n_nt* : number of nucleotides

*n_aa* : number of amino acids

*n_codon* : number of codons

*nt_to_index* : dictionary mapping one-letter nucleotides to integer indices.

*index_to_nt* : inverse of *nt_to_index*.

*aa_to_index* : dictionary mapping one-letter amino acids to integer indices.

*index_to_aa* : inverse of *aa_to_index*.

*codon_to_index* : dictionary mapping non-stop codons to integer indices.

*index_to_codon* : inverse of *codon_to_index*.

*translate_by_index* : dictionary mapping codon indices to amino-acid indices.
"""

import Bio.Alphabet.IUPAC
import Bio.Seq


ALMOST_ZERO = 1.0e-5
ALMOST_ONE = 1.0 - ALMOST_ZERO
LARGE_NUMBER = 1e4

index_to_nt = dict(enumerate(sorted(Bio.Alphabet.IUPAC.IUPACUnambiguousDNA.letters)))
nt_to_index = dict([(nt, index) for (index, nt) in index_to_nt.items()])
n_nt = len(index_to_nt)
assert len(index_to_nt) == len(nt_to_index) == n_nt

index_to_aa = dict(enumerate(sorted(Bio.Alphabet.IUPAC.IUPACProtein.letters)))
aa_to_index = dict([(aa, index) for (index, aa) in index_to_aa.items()])
n_aa = len(index_to_aa)
assert len(index_to_aa) == len(aa_to_index) == n_aa

codon_to_index = {}
index_to_codon = {}
translate_by_index = {}
i = 0
for nt1 in nt_to_index.keys():
    for nt2 in nt_to_index.keys():
        for nt3 in nt_to_index.keys():
            codon = nt1 + nt2 + nt3
            aa = str(Bio.Seq.Seq(codon).translate())
            if aa != '*':
                codon_to_index[codon] = i
                index_to_codon[i] = codon
                translate_by_index[i] = aa_to_index[aa]
                i += 1
n_codon = len(codon_to_index)
assert len(codon_to_index) == len(index_to_codon) == len(translate_by_index) == n_codon
