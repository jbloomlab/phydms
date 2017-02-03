"""Tests branch scaling.

Makes sure we can correctly re-scale branch lengths into
units of substitutions per site.

Written by Jesse Bloom.
"""

import os
import sys
import scipy
import math
import unittest
import random
import io
import copy
import phydmslib.models
import phydmslib.treelikelihood
from phydmslib.constants import *
import Bio.SeqIO
import Bio.Phylo
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


class test_branchScale_ExpCM(unittest.TestCase):
    """Tests `branchScale` of `ExpCM_empirical_phi` model."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi

    def test_branchScale(self):
        """Simulate evolution, ensure scaled branches match number of subs."""

        scipy.random.seed(1)
        random.seed(1)

        # define model, only free parameter is mu for testing simulations
        nsites = 200
        prefs = []
        minpref = 0.01
        for r in range(nsites):
            rprefs = scipy.random.dirichlet([1] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
        kappa = 4.2
        omega = 0.4
        beta = 1.5
        mu = 0.3
        if self.MODEL == phydmslib.models.ExpCM:
            phi = scipy.random.dirichlet([7] * N_NT)
            model = phydmslib.models.ExpCM(prefs, kappa=kappa, omega=omega,
                    beta=beta, mu=mu, phi=phi, freeparams=['mu'])
            partitions = pyvolvePartitionsExpCM(model)
        elif self.MODEL == phydmslib.models.ExpCM_empirical_phi:
            g = scipy.random.dirichlet([7] * N_NT)
            model = phydmslib.models.ExpCM_empirical_phi(prefs, g,
                    kappa=kappa, omega=omega, beta=beta, mu=mu,
                    freeparams=['mu'])
            partitions = pyvolvePartitionsExpCM(model)
        elif self.MODEL == phydmslib.models.YNGKP_M0:
            e_pw = scipy.asarray([scipy.random.dirichlet([7] * N_NT) for i 
                    in range(3)])
            model = phydmslib.models.YNGKP_M0(e_pw, nsites)
            partitions = pyvolvePartitionsYNGKP_M0(model)
        else:
            raise ValueError("Invalid MODEL: {0}".format(type(self.MODEL)))

        # tree is two sequences separated by a single branch
        t = 0.04 / model.branchScale
        newicktree = '(tip1:{0},tip2:{0});'.format(t / 2.0)
        pyvolvetree = pyvolve.read_tree(tree=newicktree)
        temptree = '_temp.tree'
        with open(temptree, 'w') as f:
            f.write(newicktree)
        biotree = Bio.Phylo.read(temptree, 'newick')
        os.remove(temptree)

        # Simulate evolution of two sequences separated by a long branch.
        # Then estimate subs per site in a heuristic way that will be 
        # roughly correct for short branches. Do this all several times
        # and average results to get better accuracy.
        alignment = '_temp_simulatedalignment.fasta'
        info = '_temp_info.txt'
        rates = '_temp_ratefile.txt'
        evolver = pyvolve.Evolver(partitions=partitions, tree=pyvolvetree)
        nsubs = 0 # subs in simulated seqs (estimate from Hamming distance)
        treedist = 0.0 # distance inferred by `TreeLikelihood`
        nreplicates = 300
        for i in range(nreplicates):
            evolver(seqfile=alignment, infofile=info, ratefile=rates)
            a = [(s.description, str(s.seq)) for s in Bio.SeqIO.parse(
                    alignment, 'fasta')]
            assert len(a[0][1]) == len(a[1][1]) == nsites * 3
            for f in [alignment, info, rates]:
                os.remove(f)
            for r in range(nsites):
                codon1 = a[0][1][3 * r : 3 * r + 3]
                codon2 = a[1][1][3 * r : 3 * r + 3]
                nsubs += len([j for j in range(3) if codon1[j] != codon2[j]])
            tl = phydmslib.treelikelihood.TreeLikelihood(biotree, a, 
                    copy.deepcopy(model)) 
            tl.maximizeLikelihood()
            treedist += sum([n.branch_length for n in tl.tree.get_terminals()])
        nsubs /= float(nsites * nreplicates)
        treedist /= float(nreplicates)

        # We expect nsubs = branchScale * t, but build in some tolerance
        # of rtol = 0.1 since we simulated finite number of sites.
        self.assertTrue(scipy.allclose(nsubs, model.branchScale * t, rtol=0.1),
                ("Simulated subs per site of {0} is not close "
                "to expected value of {1} (branchScale = {2}, t = {3})").format(
                nsubs, t * model.branchScale, model.branchScale, t))
        self.assertTrue(scipy.allclose(treedist, nsubs, rtol=0.1), (
                "Simulated subs per site of {0} is not close to inferred "
                "branch length of {1}").format(nsubs, treedist))

    def tearDown(self):
        """Remove some files made by `pyvolve`."""
        for f in ['custom_matrix_frequencies.txt']:
            if os.path.isfile(f):
                os.remove(f)


class test_branchScale_YNGKP_M0(test_branchScale_ExpCM):
    """Tests `branchScale` of `YNGKP_M0` model."""
    MODEL = phydmslib.models.YNGKP_M0


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
