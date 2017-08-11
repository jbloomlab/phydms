"""Functions for performing simulations, mostly using ``pyvolve``.

Written by Jesse Bloom and Sarah Hilton.
"""

import os
import sys
import scipy
import math
import phydmslib.models
from phydmslib.constants import *
import pyvolve
from tempfile import mkstemp
import random


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
                        if type(model) == phydmslib.models.ExpCM_empirical_phi_divpressure:
                            fxy *= model.omega * (1 + model.omega2 * model.deltar[r])
                        elif r + 1 in divsites:
                            fxy *= divomega
                        else:
                            fxy *= model.omega
                    if type(model) in [phydmslib.models.ExpCM,
                            phydmslib.models.ExpCM_empirical_phi, phydmslib.models.ExpCM_empirical_phi_divpressure]:
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

def simulateAlignment(model, treeFile, alignmentPrefix, randomSeed=False):
    """
    Simulate an alignment given a model and tree (units = subs/site).

    Simulations done using `pyvolve`.

    Args:
        `model` (`phydmslib.models.Models` object)
            The model used for the simulations. Only
            models that can be passed to `pyvolve.Partitions`
            are supported.
        `treeFile` (str)
            Name of newick file used to simulate the sequences.
            The branch lengths should be in substitutions per site,
            which is the default units for all `phydms` outputs.
        `alignmentPrefix`
            Prefix for the files created by `pyvolve`.

    The result of this function is a simulated FASTA alignment
    file with the name having the prefix giving by `alignmentPrefix`
    and the suffix `'_simulatedalignment.fasta'`.
    """
    if randomSeed == False:
        pass
    else:
        random.seed(randomSeed)
        
    #Transform the branch lengths by dividing by the model `branchScale`
    tree = Bio.Phylo.read(treeFile, 'newick')
    for node in tree.get_terminals() + tree.get_nonterminals():
        if (node.branch_length == None) and (node == tree.root):
            node.branch_length = 1e-06
        else:
            node.branch_length /= model.branchScale
    fd, temp_path = mkstemp()
    Bio.Phylo.write(tree, temp_path, 'newick')
    os.close(fd)
    pyvolve_tree = pyvolve.read_tree(file=temp_path)
    os.remove(temp_path)


    #Make the `pyvolve` partition
    partitions = pyvolvePartitions(model)

    #Simulate the alignment
    alignment = '{0}_simulatedalignment.fasta'.format(alignmentPrefix)
    info = '_temp_{0}info.txt'.format(alignmentPrefix)
    rates = '_temp_{0}_ratefile.txt'.format(alignmentPrefix)
    evolver = pyvolve.Evolver(partitions=partitions, tree=pyvolve_tree)
    evolver(seqfile=alignment, infofile=info, ratefile=rates)
    for f in [rates,info, "custom_matrix_frequencies.txt"]:
        if os.path.isfile(f):
            os.remove(f)
    assert os.path.isfile(alignment)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
