"""Tests random number seeding in aligment simulation.

Makes sure the random numbering seeding gives reproducible results.

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
import glob



class test_simulateAlignmentRandomSeed_ExpCM(unittest.TestCase):
    """Tests `simulateAlignment` of `simulate.py` module with different seed numbers."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi

    def setUp(self):
        """Set up parameters for test."""
        # define model
        self.nsites = 100
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

    def test_pyvovle_randomSeed(self):
        """Simulate evolution, ensure scaled branches match number of subs."""
        counter = 0
        seed = 1
        alignments = [{}, {}, {}]
        # alignments with the same seed number should be the same
        # make two alignments with the same seed number
        for counter in range(2):
            alignmentPrefix = "test_counter{0}_seed{1}".format(counter,seed)
            phydmslib.simulate.simulateAlignment(self.model, self.tree_fname, alignmentPrefix, seed)
            for s in Bio.SeqIO.parse("test_counter{0}_seed{1}_simulatedalignment.fasta".format(counter,seed), "fasta"):
                alignments[counter][s.id] = str(s.seq)
        # check they are the same
        for key in alignments[counter].keys():
            self.assertTrue(alignments[counter][key] == alignments[counter - 1][key])

        # alignments with different seed numbers should be different
        # make an alignment with a different seed number
        seed += 1
        counter += 1
        alignmentPrefix = "test_counter{0}_seed{1}".format(counter,seed)
        phydmslib.simulate.simulateAlignment(self.model, self.tree_fname, alignmentPrefix, seed)
        for s in Bio.SeqIO.parse("test_counter{0}_seed{1}_simulatedalignment.fasta".format(counter,seed), "fasta"):
            alignments[counter][s.id] = str(s.seq)
        # check they are different
        for key in alignments[counter].keys():
            self.assertFalse(alignments[counter][key] == alignments[counter - 1][key])

        # general clean-up
        os.remove(self.tree_fname)
        for fasta in glob.glob("test*simulatedalignment.fasta"):
            if os.path.isfile(fasta):
                os.remove(fasta)

    def test_simulator_randomSeed(self):
        """Simulate evolution with the `Simulator` class.
        Ensures scaled branches match number of subs.
        """
        seed = 1
        alignments = []
        simulator = phydmslib.simulate.Simulator(self.tree, self.model)
        for i in range(2):
            alignments.append(dict(simulator.simulate(seed)))
        alignments.append(dict(simulator.simulate(seed+1)))

        # check that the first two alignments are the same and the last
        # is different.
        for key in alignments[0].keys():
            self.assertTrue(alignments[0][key] == alignments[1][key])
            self.assertFalse(alignments[0][key] == alignments[2][key])

    def test_simulator_randomSeed_outside(self):
        """Simulate evolution with the `Simulator` class.
        Sets seed before calling `simulate`.
        Ensures scaled branches match number of subs.
        """
        seed = 1
        alignments = []
        simulator = phydmslib.simulate.Simulator(self.tree, self.model)

        scipy.random.seed(seed)
        for i in range(2):
            alignments.append(dict(simulator.simulate()))

        scipy.random.seed(seed)
        for i in range(2):
            alignments.append(dict(simulator.simulate()))

        # check that the first/third and second/fourth are the same
        # check that the first and second are different
        for key in alignments[0].keys():
            self.assertTrue(alignments[0][key] == alignments[2][key])
            self.assertTrue(alignments[1][key] == alignments[3][key])
            self.assertFalse(alignments[0][key] == alignments[1][key])


class test_simulateAlignmentRandomSeed_YNGKP_M0(test_simulateAlignmentRandomSeed_ExpCM):
    """Tests `simulateAlignment` of `YNGKP_M0` model."""
    MODEL = phydmslib.models.YNGKP_M0


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
