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
    DISTRIBUTIONMODEL = None

    def testExpCM_spielmanwr(self):
        """Test the `ExpCM` function `_spielman_wr`."""

        # create preference set
        random.seed(1)
        scipy.random.seed(1)
        nsites = 10
        prefs = []
        minpref = 0.01
        for r in range(nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))

        if self.MODEL == phydmslib.models.ExpCM:
            paramvalues = {
                    'eta':scipy.random.dirichlet([5] * (N_NT - 1)),
                    'omega':0.7,
                    'kappa':2.5,
                    'beta':1.2,
                    'mu':0.5,
                    }
            model = self.MODEL(prefs)
            assert set(paramvalues.keys()) == set(model.freeparams), (
                    "{0} vs {1}".format(set(paramvalues.keys()),
                    set(model.freeparams)))
            model.updateParams(paramvalues)
        elif self.MODEL == phydmslib.models.ExpCM_empirical_phi:
            g = scipy.random.dirichlet([3] * N_NT)
            paramvalues = {
                    'omega':0.7,
                    'kappa':2.5,
                    'beta':1.2,
                    'mu':0.5,
                    }
            model = self.MODEL(prefs, g)
            assert set(paramvalues.keys()) == set(model.freeparams), (
                    "{0} vs {1}".format(set(paramvalues.keys()),
                    set(model.freeparams)))
            model.updateParams(paramvalues)
        else:
            raise ValueError("Invalid BASEMODEL: {0}".format(self.MODEL))
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))

        if self.DISTRIBUTIONMODEL is None:
            pass
        elif (self.DISTRIBUTIONMODEL ==
                phydmslib.models.GammaDistributedOmegaModel):
            self.MODEL = self.DISTRIBUTIONMODEL(model, ncats=4)
        elif (self.DISTRIBUTIONMODEL ==
                phydmslib.models.GammaDistributedBetaModel):
            self.MODEL = self.DISTRIBUTIONMODEL(model, ncats=4)
        else:
            raise ValueError("Invalid DISTRIBUTIONMODEL: {0}".format(self.DISTRIBUTIONMODEL))

    def check_speilmanwr(self):
        # test `_spielman_wr` calculation
        wr = []
        for n in range(self.model.nsites):
            numerator = 0
            denominator = 0
            for x in range(N_CODON):
                for y in range(N_CODON):
                    if CODON_SINGLEMUT[x][y] and CODON_NONSYN[x][y]:
                        prx = self.model.stationarystate[n][x]
                        Prxy = self.model.Prxy[n][x][y]
                        Qxy = self.model.Qxy[x][y]
                        numerator += prx * Prxy
                        denominator += prx * Qxy
            wr.append(numerator/denominator)
        wr = scipy.array(wr)
        self.assertTrue(scipy.allclose(wr, scipy.array(self.model._spielman_wr()), rtol=0.01))

class test_empirical_phi_spielmanwr(testExpCM_spielmanwr):
    """Test the calculation of `spielmanwr` using the model `ExpCM_empirical_phi`"""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi
    DISTRIBUTIONMODEL = None

class test_empirical_phi_gammaomega_spielmanwr(testExpCM_spielmanwr):
    """Test the calculation of `spielmanwr` using the model `ExpCM_empirical_phi_gammaomega`"""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedOmegaModel

class test_expcm_gammaomega_spielmanwr(testExpCM_spielmanwr):
    """Test the calculation of `spielmanwr` using the model `ExpCM_gammaomega`"""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedOmegaModel

class test_expcm_gammabeta_spielmanwr(testExpCM_spielmanwr):
    """Test the calculation of `spielmanwr` using the model `ExpCM_gammabeta`"""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedBetaModel

class test_empirical_phi_gammabeta_spielmanwr(testExpCM_spielmanwr):
    """Test the calculation of `spielmanwr` using the model `ExpCM_empirical_phi_gammabeta`"""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    MODEL = phydmslib.models.ExpCM_empirical_phi
    DISTRIBUTIONMODEL = phydmslib.models.GammaDistributedBetaModel


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
