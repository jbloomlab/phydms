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
    """Tests simulating an alignment using both `pyvolve` and `Simulator`"""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi

    def setUp(self):
        """Set up parameters for test."""
        # define model
        self.nsites = 1000
        prefs = []
        minpref = 0.01
        for r in range(self.nsites):
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
            self.model = phydmslib.models.ExpCM(prefs, kappa=kappa, omega=omega,
                                                beta=beta, mu=mu, phi=phi,
                                                freeparams=['mu'])
        elif self.MODEL == phydmslib.models.ExpCM_empirical_phi:
            g = scipy.random.dirichlet([7] * N_NT)
            self.model = phydmslib.models.ExpCM_empirical_phi(prefs, g,
                                                              kappa=kappa,
                                                              omega=omega,
                                                              beta=beta,
                                                              mu=mu,
                                                              freeparams=['mu'])
        elif self.MODEL == phydmslib.models.YNGKP_M0:
            e_pw = scipy.asarray([scipy.random.dirichlet([7] * N_NT) for i
                                  in range(3)])
            self.model = phydmslib.models.YNGKP_M0(e_pw, self.nsites)
        else:
            raise ValueError("Invalid MODEL: {0}".format(type(self.MODEL)))

        # make a test tree
        # tree is two sequences separated by a single branch
        # the units are in sub/site
        self.t = 0.04
        newicktree = '(tip1:{0},tip2:{0});'.format(self.t / 2.0)
        self.tree_fname = '_temp.tree'
        with open(self.tree_fname, 'w') as f:
            f.write(newicktree)
        self.tree = Bio.Phylo.read(self.tree_fname, 'newick')

        self.rescaled_tree = copy.deepcopy(self.tree)
        # re-scale the branch lengths
        for node in self.rescaled_tree.get_terminals() + self.rescaled_tree.get_nonterminals():
            if node.branch_length:
                node.branch_length /= self.model.branchScale

    def test_simulateAlignment_Simulator(self):
        """Simulate evolution with the `Simulator` class.
        Ensures scaled branches match number of subs.
        """
        # set random seeds
        scipy.random.seed(1)
        random.seed(1)
        simulator = phydmslib.simulate.Simulator(self.tree, self.model)
        alignment = simulator.simulate(1)
        self.check_alignment(alignment)

    def test_simulateAlignment_pyvolve(self):
        """Simulate evolution with `pyvolve` class.
        Ensures scaled branches match number of subs.
        """
        # set random seeds
        scipy.random.seed(1)
        random.seed(1)

        # simulate the alignment using `pyvlove`
        alignment = phydmslib.simulate.simulateAlignment(self.model,
                                                         self.tree_fname,
                                                         "{0}_pyvolve".format("test"),
                                                         randomSeed=1)
        os.remove(self.tree_fname)
        assert len(alignment) == 1, "Only expected one `pyvolve` alignment"

        a = [(s.description, str(s.seq)) for s in Bio.SeqIO.parse(
                alignment[0], 'fasta')]
        self.check_alignment(a)
        if os.path.isfile(alignment[0]):
            os.remove(alignment[0])

    def check_alignment(self, a):
        # check and see if the simulated alignment has the expected number of
        # subs exists
        nsubs = 0  # subs in simulated seqs (estimate from Hamming distance)
        treedist = 0.0  # distance inferred by `TreeLikelihood`
        assert len(a[0][1]) == len(a[1][1]) == self.nsites * 3
        for r in range(self.nsites):
            codon1 = a[0][1][3 * r: 3 * r + 3]
            codon2 = a[1][1][3 * r: 3 * r + 3]
            nsubs += len([j for j in range(3) if codon1[j] != codon2[j]])
        nsubs /= float(self.nsites)
        tl = phydmslib.treelikelihood.TreeLikelihood(self.rescaled_tree,
                                                     a, self.model)
        tl.maximizeLikelihood()
        treedist += sum([n.branch_length for n in tl.tree.get_terminals()])

        # We expect nsubs = t, but build in some tolerance
        # with rtol since we simulated finite number of sites.
        self.assertTrue(scipy.allclose(nsubs, self.t, rtol=0.2),
                        ("Simulated subs per site of {0} is not close "
                         "to expected value of {1} (branchScale = {2},"
                         " t = {3})").format(
                nsubs, self.t, self.model.branchScale, self.t))
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
