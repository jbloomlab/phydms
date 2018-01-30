"""Tests the calculation of spielmanwr, following Spielman and Wilke, 2015.

Written by Jesse Bloom and Sarah Hilton."""


import random
import unittest
import scipy
from phydmslib.constants import *
import phydmslib.models

class testExpCM_spielmanwr(unittest.TestCase):
    """Test the calculation of `spielmanwr` using the model `ExpCM`."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM

    def testExpCM_spielmanwr(self):
        """Test the `ExpCM` function `_spielman_wr`."""

        # create models
        random.seed(1)
        scipy.random.seed(1)
        nsites = 10
        g = scipy.random.dirichlet([5] * N_NT)
        prefs = []
        minpref = 0.01
        for r in range(nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))

        if self.MODEL == phydmslib.models.ExpCM:
            self.model = phydmslib.models.ExpCM(prefs)
        elif self.MODEL == phydmslib.models.ExpCM_empirical_phi:
            self.model = phydmslib.models.ExpCM_empirical_phi(prefs, g)
        else:
            raise ValueError("Invalid MODEL: {0}".format(self.MODEL))

        # test `_spielman_wr` calculation
        wr = []
        for n in range(self.model.nsites):
            numerator = 0
            denominator = 0
            for x in range(N_CODON):
                for y in range(N_CODON):
                    if CODON_SINGLEMUT[x][y] and CODON_NONSYN[x][y]:
                        prx = self.model.stationarystate[n][x]
                        Prxy = (self.model.Prxy[n][x][y])
                        Qxy = self.model.Qxy[x][y]
                        numerator += prx * Prxy
                        denominator += prx * Qxy
            wr.append(numerator/denominator)
        wr = scipy.array(wr)
        self.assertTrue(scipy.allclose(wr, scipy.array(self.model.spielman_wr(norm=False)), rtol=0.01))
        self.assertTrue(scipy.allclose(wr/self.model.omega, scipy.array(self.model.spielman_wr()), rtol=0.01))


class test_empirical_phi_spielmanwr(testExpCM_spielmanwr):
    """Test the calculation of `spielmanwr` using the model `ExpCM_empirical_phi`"""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
