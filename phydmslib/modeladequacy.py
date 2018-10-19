"""Functions for performing model adequacy tests.

Written by Sarah Hilton.
"""

from phydmslib.constants import *
import phydmslib.models
import pandas as pd
import scipy
import math
import pandas as pd

scipy.random.seed(0)

def translate_with_gaps(seq):
    """Translate from nucleotides to amino acids.

    Args:
        `seq` (str)
            The nucleotide sequence.
    Returns:
        The amino-acid sequence

    >>> s1 = "ATGATG"
    >>> s2 = "CTT---ATG"
    >>> translate_with_gaps(s1) == "MM"
    True
    >>> translate_with_gaps(s2) == "L-M"
    True

    """
    assert len(seq) % 3 == 0, "Sequence is not divisible by 3."
    prot_seq = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if codon == "---":
            AA = "-"
        else:
            codon = CODON_TO_INDEX[codon]
            AA = INDEX_TO_AA[CODON_TO_AA[codon]]
        prot_seq.append(AA)
    return "".join(prot_seq)


def calc_aa_frequencies(alignment):
    """Calculate amino-acid frequencies from a codon alignment.

    Args:
        `alignment` (list)
            Alignment of codon sequences as a list of tuples, (seq_id, seq)
    Returns:
        `pandas` dataframe of amino-acid frequencies by site

    >>> answer = pd.DataFrame({"site": [1, 2], "A": [0.0, 0.0],\
                               "C": [0.0, 0.0], "D": [0.0, 0.0],\
                               "E": [0.0, 0.0], "F": [0.0, 0.0],\
                               "G": [0.0, 0.0], "H": [0.0, 0.0],\
                               "I": [0.0, 0.0], "K": [0.0, 0.0],\
                               "L": [0.0, 0.0], "M": [0.0, 0.0],\
                               "N": [0.0, 0.0], "P": [0.0, 0.0],\
                               "Q": [0.0, 0.0], "R": [0.0, 0.0],\
                               "S": [0.0, 0.0], "T": [0.0, 0.0],\
                               "V": [0.0, 0.0], "W": [0.0, 0.0],\
                               "Y": [0.0, 0.0]}, columns=["site","A","C",\
                                                          "D","E","F","G",\
                                                          "H","I","K","L",\
                                                          "M","N","P","Q",\
                                                          "R","S","T","V",\
                                                          "W","Y"])
    >>> align1 = [("seq_1", "ATGATG"), ("seq_2", "CTTATG")]
    >>> align2 = [("seq_1", "ATGATG"), ("seq_2", "CTT---")]
    >>> answer1 = answer.copy()
    >>> answer1[["L", "M"]] = pd.DataFrame({"L": [0.5, 0.0],\
                                            "M": [0.5, 1.0]})
    >>> answer1.equals(calc_aa_frequencies(align1))
    True
    >>> answer1.equals(calc_aa_frequencies(align2))
    True

    """
    # Read in the alignnment
    assert scipy.all(scipy.array([len(s[1]) % 3 for s in alignment]) == 0),\
        "At least one sequence in the alignment is not a multiple of 3."
    seqlength = len(alignment[0][1]) // 3
    df = {k: [0 for x in range(seqlength)] for k in list(INDEX_TO_AA.values())}

    # count amino acid frequencies
    for seq in alignment:
        for i, aa in enumerate(translate_with_gaps(seq[1])):
            if aa != '-':
                df[aa][i] += 1
    df = pd.DataFrame(df)

    # # Normalize the dataframe
    df = df.div(df.sum(axis=1), axis=0)
    assert scipy.allclose(df.sum(axis=1), 1, atol=0.005)

    # final formatting
    aa = [x for x in INDEX_TO_AA.values()]
    aa.sort()  # ABC order
    final_cols = ["site"]
    final_cols.extend(aa)
    df["site"] = [x+1 for x in range(len(df))]
    df = df[final_cols]
    return df


def prefDistance(pi1, pi2, distmetric):
    """Compute the distance between two arrays of preferences.

    Args:
        `pi1` and `pi2` (array-like)
            Two arrays of preferences.
        `distmetric` (string)
            Distance metric to use. Can be:
                - `half_sum_abs_diff`: half sum absolute value of difference
                - `JensenShannon`: square root of Jensen-Shannon divergence

    Returns:
        The distance between `pi1` and `pi2`.

    >>> pi1 = [0.5, 0.2, 0.3]
    >>> pi2 = [0.2, 0.4, 0.4]
    >>> scipy.allclose(prefDistance(pi1, pi1, 'half_sum_abs_diff'), 0)
    True
    >>> scipy.allclose(prefDistance(pi1, pi1, 'JensenShannon'), 0)
    True
    >>> scipy.allclose(prefDistance(pi1, pi2, 'half_sum_abs_diff'), 0.3)
    True
    >>> scipy.allclose(prefDistance(pi1, pi2, 'JensenShannon'), 0.2785483)
    True

    """
    pi1 = scipy.array(pi1)
    pi2 = scipy.array(pi2)
    assert len(pi1) == len(pi2)
    assert scipy.allclose(pi1.sum(), 1, atol=0.005)
    assert scipy.allclose(pi2.sum(), 1, atol=0.005)
    assert scipy.all(pi1 >= 0)
    assert scipy.all(pi2 >= 0)

    if distmetric == 'half_sum_abs_diff':
        dist = (scipy.fabs(pi1 - pi2)).sum() / 2.0

    elif distmetric == 'JensenShannon':
        dist = math.sqrt(divJensenShannon(pi1, pi2))

    else:
        raise ValueError('Invalid `distmetric` {0}'.format(distmetric))

    return dist


def divJensenShannon(p1, p2):
    """Jensen-Shannon divergence between two distributions.

    The logarithms are taken to base 2, so the result will be
    between 0 and 1.
    Args:
        `p1` and `p2` (array-like)
            The two distributions for which we compute divergence.
    Returns:
        The Jensen-Shannon divergence as a float.
    >>> p1 = [0.5, 0.2, 0.2, 0.1]
    >>> p2 = [0.4, 0.1, 0.3, 0.2]
    >>> p3 = [0.0, 0.2, 0.2, 0.6]
    >>> scipy.allclose(divJensenShannon(p1, p1), 0, atol=1e-5)
    True
    >>> scipy.allclose(divJensenShannon(p1, p2), 0.035789, atol=1e-5)
    True
    >>> scipy.allclose(divJensenShannon(p1, p3), 0.392914, atol=1e-5)
    True

    """
    p1 = scipy.array(p1)
    p2 = scipy.array(p2)

    def _kldiv(a, b):
        with scipy.errstate(all='ignore'):
            kl = scipy.sum([v for v in a * scipy.log2(a / b)
                            if not scipy.isnan(v)])
        return kl

    m = 0.5 * (p1 + p2)

    return 0.5 * (_kldiv(p1, m) + _kldiv(p2, m))

# Note from SKH to JDB
# I don't use this function in my script but it seemed like it might be useful
# at some point in someone's analyis
def make_stationary_state_prefs(model):
    """Create a preference set from stationary state amino-acid frequencies.

    Args:
        `model` (phydmslib.models.ExpCM)

    Returns:
        `prefs` (`pandas` dataframe)
            The stationary state amino-acid frequencies in preference format.

    """
    header = [INDEX_TO_AA[x] for x in range(20)]
    frequencies = calc_stationary_state_freqs(model)
    prefs = pd.DataFrame(frequencies, columns=header)
    prefs.insert(0, "site", [x+1 for x in range(len(prefs))])
    return prefs


def calc_stationary_state_freqs(model):
    """Calculate the stationary state amino-acids frequencies of a model.

    Args:
        `model` (phydmslib.models.ExpCM)

    Returns:
        frequencies (`numpy.ndarray` of floats)
            The stationary state amino-acid frequencies.
            frequencies[r][a] is the statioanry state frequence of amino acid
            `a` at site `r`.

    """
    def _calc_aa_freq(aa, ss):
        """Calc the frequency of a single site/amino acid pair."""
        codon_indices = scipy.where(CODON_TO_AA == aa)[0]
        return ss[codon_indices].sum()

    frequencies = []
    for r in range(model.nsites):
        aminoacid_ss = scipy.array([_calc_aa_freq(aa, model.stationarystate[r])
                                    for aa in range(N_AA)])
        frequencies.append(aminoacid_ss)
    return scipy.array(frequencies)


def make_expcm(model_fname, prefs):
    """Make an ExpCM from a model params file."""
    params = pd.read_csv(model_fname, engine="python", sep=" = ", header=None)
    params = dict(zip(params[0], params[1]))
    params["phiT"] = 1 - sum([params[x] for x in params.keys()
                             if x.startswith("phi")])
    phi = scipy.array([params["phi{0}".format(INDEX_TO_NT[x])]
                       for x in range(N_NT)])
    return phydmslib.models.ExpCM(prefs, kappa=params["kappa"],
                                  omega=params["omega"], beta=params["beta"],
                                  mu=0.3, phi=phi, freeparams=['mu'])


if __name__ == '__main__':
    import doctest
    doctest.testmod()
