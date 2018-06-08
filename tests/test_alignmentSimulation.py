"""Tests alignment simulation.

Makes sure we can correctly simulate an alignment given a model and a tree.

Written by Sarah Hilton and Jesse Bloom.
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
import phydmslib.simulate
from phydmslib.constants import *
import Bio.SeqIO
import Bio.Phylo
import pyvolve



class test_simulateAlignment_ExpCM(unittest.TestCase):
    """Tests `simulateAlignment` of `simulate.py` module."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi

    def test_simulateAlignment(self):
        """Simulate evolution, ensure scaled branches match number of subs."""

        scipy.random.seed(1)
        random.seed(1)

        alignmentPrefix = "test"

        # define model
        nsites = 1000
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
        elif self.MODEL == phydmslib.models.ExpCM_empirical_phi:
            g = scipy.random.dirichlet([7] * N_NT)
            model = phydmslib.models.ExpCM_empirical_phi(prefs, g,
                    kappa=kappa, omega=omega, beta=beta, mu=mu,
                    freeparams=['mu'])
        elif self.MODEL == phydmslib.models.YNGKP_M0:
            e_pw = scipy.asarray([scipy.random.dirichlet([7] * N_NT) for i
                    in range(3)])
            model = phydmslib.models.YNGKP_M0(e_pw, nsites)
        else:
            raise ValueError("Invalid MODEL: {0}".format(type(self.MODEL)))

        # make a test tree
        # tree is two sequences separated by a single branch
        # the units are in sub/site
        t = 0.04
        newicktree = '(tip1:{0},tip2:{0});'.format(t / 2.0)
        temptree = '_temp.tree'
        with open(temptree, 'w') as f:
            f.write(newicktree)

        # simulate the alignment
        phydmslib.simulate.simulateAlignment(model, temptree, alignmentPrefix)

        # read in the test tree, re-scale the branch lengths, and remove the file
        biotree = Bio.Phylo.read(temptree, 'newick')
        os.remove(temptree)
        for node in biotree.get_terminals() + biotree.get_nonterminals():
            if node.branch_length:
                node.branch_length /= model.branchScale

        # check and see if the simulated alignment has the expected number of
        # subs exists
        alignment = '{0}_simulatedalignment.fasta'.format(alignmentPrefix)
        nsubs = 0 # subs in simulated seqs (estimate from Hamming distance)
        treedist = 0.0 # distance inferred by `TreeLikelihood`
        a = [(s.description, str(s.seq)) for s in Bio.SeqIO.parse(
                alignment, 'fasta')]
        assert len(a[0][1]) == len(a[1][1]) == nsites * 3
        for f in [alignment]:
            if os.path.isfile(f):
                os.remove(f)
        for r in range(nsites):
            codon1 = a[0][1][3 * r : 3 * r + 3]
            codon2 = a[1][1][3 * r : 3 * r + 3]
            nsubs += len([j for j in range(3) if codon1[j] != codon2[j]])
        nsubs /= float(nsites)
        tl = phydmslib.treelikelihood.TreeLikelihood(biotree, a, model)
        tl.maximizeLikelihood()
        treedist += sum([n.branch_length for n in tl.tree.get_terminals()])

        # We expect nsubs = t, but build in some tolerance
        # with rtol since we simulated finite number of sites.
        self.assertTrue(scipy.allclose(nsubs, t, rtol=0.2),
                ("Simulated subs per site of {0} is not close "
                "to expected value of {1} (branchScale = {2}, t = {3})").format(
                nsubs, t, model.branchScale, t))
        self.assertTrue(scipy.allclose(treedist, nsubs, rtol=0.2), (
                "Simulated subs per site of {0} is not close to inferred "
                "branch length of {1}").format(nsubs, treedist))

    def tearDown(self):
        """Remove some files made by `pyvolve`."""
        for f in ['custom_matrix_frequencies.txt']:
            if os.path.isfile(f):
                os.remove(f)


class test_simulateAlignment_YNGKP_M0(test_simulateAlignment_ExpCM):
    """Tests `simulateAlignment` of `YNGKP_M0` model."""
    MODEL = phydmslib.models.YNGKP_M0


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
