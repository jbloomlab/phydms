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


def pyvolvePartitions(model, divselection=None):
    """Get list of `pyvolve` partitions for `model`.

    Args:
        `model` (`phydmslib.models.Models` object)
            The model used for the simulations. Currently only
            certain `Models` are supported (e.g., `YNGKP`,
            `ExpCM`)
        `divselection` (`None` or 2-tuple `(divomega, divsites)`)
            Set this option if you want to simulate a subset of sites
            as under diversifying selection (e.g., an `omega` different
            than that used by `model`. In this case, `divomega` is
            the omega for this subset of sites, and `divsites` is a list
            of the sites in 1, 2, ... numbering.

    Returns:
        `partitions` (`list` of `pyvolve.Partition` objects)
            Can be fed into `pyvolve.Evolver` to simulate evolution.
    """
    codons = pyvolve.genetics.Genetics().codons
    codon_dict = pyvolve.genetics.Genetics().codon_dict
    pyrims = pyvolve.genetics.Genetics().pyrims
    purines = pyvolve.genetics.Genetics().purines

    if divselection:
        (divomega, divsites) = divselection
    else:
        divsites = []

    assert all([1 <= r <= model.nsites for r in divsites])

    partitions = []
    for r in range(model.nsites):
        matrix = scipy.zeros((len(codons), len(codons)), dtype='float')
        for (xi, x) in enumerate(codons):
            for (yi, y) in enumerate(codons):
                ntdiffs = [(x[j], y[j]) for j in range(3) if x[j] != y[j]]
                if len(ntdiffs) == 1:
                    (xnt, ynt) = ntdiffs[0]
                    qxy = 1.0
                    if (xnt in purines) == (ynt in purines):
                        qxy *= model.kappa
                    (xaa, yaa) = (codon_dict[x], codon_dict[y])
                    fxy = 1.0
                    if xaa != yaa:
                        if r + 1 in divsites:
                            fxy *= divomega
                        else:
                            fxy *= model.omega
                    if type(model) in [phydmslib.models.ExpCM,
                            phydmslib.models.ExpCM_empirical_phi]:
                        qxy *= model.phi[NT_TO_INDEX[ynt]]
                        pix = model.pi[r][AA_TO_INDEX[xaa]]**model.beta
                        piy = model.pi[r][AA_TO_INDEX[yaa]]**model.beta
                        if abs(pix - piy) > ALMOST_ZERO:
                            fxy *= math.log(piy / pix) / (1.0 - pix / piy)
                    elif type(model) == phydmslib.models.YNGKP_M0:
                        for p in range(3):
                            qxy *= model.phi[p][NT_TO_INDEX[y[p]]]
                    else:
                        raise ValueError("Can't handle model type {0}".format(
                                type(model)))
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
