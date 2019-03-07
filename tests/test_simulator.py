"""Tests `Simulator` class to ensure reproducible results.

Compares an alignment created by `Simulator` to expected results.

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
from phydmslib.file_io import ReadCodonAlignment
import Bio.SeqIO
import Bio.Phylo
import pyvolve


class test_Simulator_ExpCM(unittest.TestCase):
    """Tests the `Simulator` class`"""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi
    EXPECTED = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            'expected_simulator_results/expected_simulator_ExpCM.fasta')

    def setUp(self):
        """Set up parameters for test."""
        # define model
        scipy.random.seed(0)
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
        elif self.MODEL == phydmslib.models.ExpCM_empirical_phi_divpressure:
            divpressure = scipy.ndarray(1000, dtype=float)
            for x in range(0, 500):
                divpressure[x] = 1
            scipy.random.shuffle(divpressure)
            g = scipy.random.dirichlet([7] * N_NT)
            self.model = phydmslib.models.ExpCM_empirical_phi_divpressure(
                prefs, g, kappa=kappa, omega=omega, omega2=self.OMEGA_2,
                beta=beta, mu=mu, divPressureValues=divpressure, freeparams=['mu'])
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

    def test_simulateAlignment_Simulator(self):
        """Compare simulation to known sequence."""
        # simulate an alignment with seeds 1 and 2
        simulator = phydmslib.simulate.Simulator(self.tree, self.model)
        seed1 = dict(simulator.simulate(1))
        seed2 = dict(simulator.simulate(2))

        # read in expected results
        expected = dict(ReadCodonAlignment(self.EXPECTED,
                                           checknewickvalid=True))

        # checks
        self.assertTrue(len(seed1) == len(seed2) == len(expected),
                        "Alignments different lengths")
        for seq_id in seed1.keys():
            if self.MODEL == phydmslib.models.ExpCM_empirical_phi_divpressure:
                self.assertTrue(
                    seed1[seq_id] != expected[seq_id] and
                    seed2[seq_id] != expected[seq_id],
                    ("Sequences generated with divpressure should be"
                     "different from alignments without divpressure."))
            else:
                self.assertTrue(seed1[seq_id] == expected[seq_id],
                                ("Sequence {0} different between simulated "
                                 "and expected alignment".format(seq_id)))
                self.assertFalse(seed2[seq_id] == expected[seq_id],
                                 ("Sequence {0} from seed 2 should be "
                                  "different from expected "
                                  "alignment".format(seq_id)))


class test_Simulator_YNGKP_MO(test_Simulator_ExpCM):
    """Tests `Simulator` simulation of `YNGKP_M0` model."""
    MODEL = phydmslib.models.YNGKP_M0
    EXPECTED = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            'expected_simulator_results/expected_simulator_YNGKP_MO.fasta')

class test_Simulator_ExpCM_divpressure(test_Simulator_ExpCM):
    """Tests `Simulator` simulation of `ExpCM` model with divpressure flag."""

    MODEL = phydmslib.models.ExpCM_empirical_phi_divpressure
    OMEGA_2 = 0.25


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
