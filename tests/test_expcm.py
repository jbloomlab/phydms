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

        # make sure phi is correct
        self.assertTrue(scipy.allclose(self.phi, self.expcm.phi))

        # make sure Prxy has rows summing to zero
        self.assertFalse(scipy.isnan(self.expcm.Prxy).any())
        self.assertFalse(scipy.isinf(self.expcm.Prxy).any())
        diag = scipy.eye(N_CODON, dtype='bool')
        for r in range(self.nsites):
            self.assertTrue(scipy.allclose(0, scipy.sum(self.expcm.Prxy[r], 
                    axis=1)))
            self.assertTrue(scipy.allclose(0, self.expcm.Prxy[r].sum()))
            self.assertTrue((self.expcm.Prxy[r][diag] <= 0).all())
            self.assertTrue((self.expcm.Prxy[r][~diag] >= 0).all())

        # make sure prx sums to 1 for each r
        self.assertTrue((self.expcm.prx >= 0).all())
        for r in range(self.nsites):
            self.assertTrue(scipy.allclose(1, self.expcm.prx[r].sum()))

        # check that prx is eigenvector or Prxy for the same r,
        # but not different r
        for r in range(self.nsites):
            self.assertTrue(scipy.allclose(0, scipy.dot(self.expcm.prx[r],
                    self.expcm.Prxy[r])))
            if r > 0:
                self.assertFalse(scipy.allclose(0, scipy.dot(self.expcm.prx[r],
                        self.expcm.Prxy[r - 1])))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
