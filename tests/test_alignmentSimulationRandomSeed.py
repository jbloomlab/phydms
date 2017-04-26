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

    def test_simulateAlignmentRandomSeed(self):
        """Simulate evolution, ensure scaled branches match number of subs."""

        scipy.random.seed(1)
        random.seed(1)

        # define model
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
        t = 0.04 / model.branchScale
        newicktree = '(tip1:{0},tip2:{0});'.format(t / 2.0)
        temptree = '_temp.tree'
        with open(temptree, 'w') as f:
            f.write(newicktree)

        counter = 0
        seed = 1
        alignments = [{}, {}, {}]
        # alignments with the same seed number should be the same
        # make two alignments with the same seed number
        for counter in range(2):
            alignmentPrefix = "test_counter{0}_seed{1}".format(counter,seed)
            phydmslib.simulate.simulateAlignment(model, temptree, alignmentPrefix, seed)
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
        phydmslib.simulate.simulateAlignment(model, temptree, alignmentPrefix, seed)
        for s in Bio.SeqIO.parse("test_counter{0}_seed{1}_simulatedalignment.fasta".format(counter,seed), "fasta"):
            alignments[counter][s.id] = str(s.seq)
        # check they are different
        for key in alignments[counter].keys():
            self.assertFalse(alignments[counter][key] == alignments[counter - 1][key])


        # general clean-up
        os.remove(temptree)
        for fasta in glob.glob("test*simulatedalignment.fasta"):
            if os.path.isfile(fasta):
                os.remove(fasta)


class test_simulateAlignmentRandomSeed_YNGKP_M0(test_simulateAlignmentRandomSeed_ExpCM):
    """Tests `simulateAlignment` of `YNGKP_M0` model."""
    MODEL = phydmslib.models.YNGKP_M0


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
