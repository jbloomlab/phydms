"""Tests `phydmslib.models.ExpCM` class.

Written by Jesse Bloom."""


import unittest
import scipy
from phydmslib.constants import *
import phydmslib.models


class testExpCM(unittest.TestCase):

    def setUp(self):
        """Initialize `ExpCM` with specific values."""
        scipy.random.seed(1)
        self.omega = 0.7
        self.kappa = 2.5
        self.beta = 1.5
        self.phi = scipy.random.dirichlet([2] * N_NT)
        self.nsites = 4
        self.prefs = []
        for r in range(self.nsites):
            self.prefs.append(dict(zip(AA_TO_INDEX.keys(), 
                    scipy.random.dirichlet([0.5] * N_AA))))
        self.expcm = phydmslib.models.ExpCM(self.prefs, kappa=self.kappa,
                omega=self.omega, beta=self.beta, phi=self.phi)

    def test_initialize_ExpCM(self):
        """Make sure `ExpCM` initialized to expected values."""
        self.assertEqual(self.nsites, self.expcm.nsites)
        self.assertTrue(scipy.allclose(self.phi, self.expcm.phi))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
