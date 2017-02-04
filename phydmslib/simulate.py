"""Functions for performing simulations, mostly using ``pyvolve``.

Written by Jesse Bloom.
"""

import os
import sys
import scipy
import math
import phydmslib.models
from phydmslib.constants import *
import pyvolve


def pyvolvePartitionsYNGKP_M0(model):
    """Get list of `pyvolve` partitions for a `YNGKP_M0`.

    The resulting list of partitions can be fed directly into
    `pyvolve.Evolver` to simulation evolution of sequences
    under the `YNGKP_M0` specified by `model`.

    Args:
        `model` (`YNGKP_M0` from `phydmslib.models`)
    Returns:
        `partitions` (`list` of `pyvolve.Partition` objects)
    """
    assert type(model) == phydmslib.models.YNGKP_M0

    codons = pyvolve.genetics.Genetics().codons
    codon_dict = pyvolve.genetics.Genetics().codon_dict
    pyrims = pyvolve.genetics.Genetics().pyrims
    purines = pyvolve.genetics.Genetics().purines

    matrix = scipy.zeros((len(codons), len(codons)), dtype='float')
    for (xi, x) in enumerate(codons):
        for (yi, y) in enumerate(codons):
            ntdiffs = [(x[j], y[j]) for j in range(3) if x[j] != y[j]]
            if len(ntdiffs) == 1:
                (xnt, ynt) = ntdiffs[0]
                qxy = 1
                for p in range(3):
                    qxy *= model.phi[p][NT_TO_INDEX[y[p]]]
                assert qxy == model.Phi_x[CODON_TO_INDEX[y]]
                if (xnt in purines) == (ynt in purines):
                    qxy *= model.kappa
                (xaa, yaa) = (codon_dict[x], codon_dict[y])
                if xaa == yaa:
                    fxy = 1.0
                else:
                    fxy = model.omega
                matrix[xi][yi] = model.mu * qxy * fxy
                assert scipy.allclose(matrix[xi][yi], model.Pxy[0][CODON_TO_INDEX[x]][CODON_TO_INDEX[y]])
        matrix[xi][xi] = -matrix[xi].sum()

    # create model in way that captures annoying print statements in pyvolve
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')
    try:
        m = pyvolve.Model("custom", {"matrix":matrix})
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

    partitions = [pyvolve.Partition(models=m, size=model.nsites)]
    return partitions


def pyvolvePartitionsExpCM(model):
    """Get list of `pyvolve` partitions for an `ExpCM`.

    The resulting list of partitions can be fed directly into 
    `pyvolve.Evolver` to simulate evolution of sequences under
    the `ExpCM` specified by `model`.

    Args:
        `model` (`ExpCM` or `ExpCM_empirical_phi` from `phydmslib.models`)
    Returns:
        `partitions` (`list` of `pyvolve.Partition` objects)
    """
    assert type(model) in [phydmslib.models.ExpCM, 
            phydmslib.models.ExpCM_empirical_phi]

    codons = pyvolve.genetics.Genetics().codons
    codon_dict = pyvolve.genetics.Genetics().codon_dict
    pyrims = pyvolve.genetics.Genetics().pyrims
    purines = pyvolve.genetics.Genetics().purines

    partitions = []
    for r in range(model.nsites):
        matrix = scipy.zeros((len(codons), len(codons)), dtype='float')
        for (xi, x) in enumerate(codons):
            for (yi, y) in enumerate(codons):
                ntdiffs = [(x[j], y[j]) for j in range(3) if x[j] != y[j]]
                if len(ntdiffs) == 1:
                    (xnt, ynt) = ntdiffs[0]
                    if (xnt in purines) == (ynt in purines):
                        qxy = model.kappa * model.phi[NT_TO_INDEX[ynt]]
                    else:
                        qxy = model.phi[NT_TO_INDEX[ynt]]
                    (xaa, yaa) = (codon_dict[x], codon_dict[y])
                    if xaa == yaa:
                        fxy = 1.0
                    else:
                        pix = model.pi[r][AA_TO_INDEX[xaa]]**model.beta
                        piy = model.pi[r][AA_TO_INDEX[yaa]]**model.beta
                        if abs(pix - piy) < 1e-6:
                            fxy = model.omega
                        else:
                            fxy = model.omega * math.log(piy / pix) / (1.0 - pix / piy)
                    matrix[xi][yi] = model.mu * qxy * fxy
            matrix[xi][xi] = -matrix[xi].sum()

        # create model in way that captures annoying print statements in pyvolve
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')
        try:
            m = pyvolve.Model("custom", {"matrix":matrix})
        finally:
            sys.stdout.close()
            sys.stdout = old_stdout

        partitions.append(pyvolve.Partition(models=m, size=1))
    return partitions



if __name__ == '__main__':
    import doctest
    doctest.testmod()
