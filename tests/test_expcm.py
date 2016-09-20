"""Tests `phydmslib.models.ExpCM` class.

Written by Jesse Bloom."""


import random
import unittest
import scipy
from phydmslib.constants import *
import phydmslib.models


class testExpCM(unittest.TestCase):

    def test_ExpCM(self):
        """Initialize `ExpCM`, test values, update, test again."""

        # create preferences
        scipy.random.seed(1)
        random.seed(1)
        self.nsites = 4
        self.prefs = []
        for r in range(self.nsites):
            self.prefs.append(dict(zip(AA_TO_INDEX.keys(), 
                    scipy.random.dirichlet([0.5] * N_AA))))

        # create initial ExpCM
        params = {'omega':0.7, 'kappa':2.5, 'beta':1.5,
                  'phi':scipy.random.dirichlet([2] * N_NT)}
        self.expcm = phydmslib.models.ExpCM(self.prefs, kappa=params['kappa'],
                omega=params['omega'], beta=params['beta'], phi=params['phi'])
        self.assertTrue(scipy.allclose(params['phi'], self.expcm.phi))
        self.check_ExpCM_attributes()

        # now update values and re-check several times
        for update in range(5):
            params = {'omega':random.uniform(ALMOST_ZERO, 5),
                      'kappa':random.uniform(1.0, 10.0),
                      'beta':random.uniform(ALMOST_ZERO, 4),
                      'eta':scipy.random.dirichlet([2] * (N_NT - 1))
                     }
            self.expcm.updateParams(omega=params['omega'], beta=params['beta'],
                    kappa=params['kappa'], eta=params['eta'])
            self.check_ExpCM_attributes()


    def check_ExpCM_attributes(self):
        """Make sure `ExpCM` has the expected attribute values."""
        self.assertEqual(self.nsites, self.expcm.nsites)

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

        # prx is eigenvector or Prxy for the same r, but not different r
        for r in range(self.nsites):
            self.assertTrue(scipy.allclose(0, scipy.dot(self.expcm.prx[r],
                    self.expcm.Prxy[r])))
            if r > 0:
                self.assertFalse(scipy.allclose(0, scipy.dot(self.expcm.prx[r],
                        self.expcm.Prxy[r - 1])))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
