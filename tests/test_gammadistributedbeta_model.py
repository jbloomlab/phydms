"""Tests `phydmslib.models.GammaDistributedBetaModel` class.

Written by Jesse Bloom and Sarah Hilton.
"""


import random
import unittest
import scipy
import scipy.linalg
import scipy.optimize
from phydmslib.constants import *
import phydmslib.models


class test_GammaDistributedBeta_ExpCM(unittest.TestCase):
    """Test gamma distributed beta for `ExpCM`."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    BASEMODEL = phydmslib.models.ExpCM

    def test_GammaDistributedBeta(self):
        """Initialize, test values, update, test again."""

        random.seed(1)
        scipy.random.seed(1)
        nsites = 10

        # create preference set
        prefs = []
        minpref = 0.01
        for r in range(nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))

        if self.BASEMODEL == phydmslib.models.ExpCM:
            paramvalues = {
                    'eta':scipy.random.dirichlet([5] * (N_NT - 1)),
                    'omega':0.7,
                    'kappa':2.5,
                    'beta':1.2,
                    'mu':0.5,
                    }
            basemodel = self.BASEMODEL(prefs)
            assert set(paramvalues.keys()) == set(basemodel.freeparams), (
                    "{0} vs {1}".format(set(paramvalues.keys()),
                    set(basemodel.freeparams)))
            basemodel.updateParams(paramvalues)
        elif self.BASEMODEL == phydmslib.models.ExpCM_empirical_phi:
            g = scipy.random.dirichlet([3] * N_NT)
            paramvalues = {
                    'omega':0.7,
                    'kappa':2.5,
                    'beta':1.2,
                    'mu':0.5,
                    }
            basemodel = self.BASEMODEL(prefs, g)
            assert set(paramvalues.keys()) == set(basemodel.freeparams), (
                    "{0} vs {1}".format(set(paramvalues.keys()),
                    set(basemodel.freeparams)))
            basemodel.updateParams(paramvalues)
        else:
            raise ValueError("Invalid BASEMODEL: {0}".format(self.BASEMODEL))
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))

        ncats = 4
        gammamodel = phydmslib.models.GammaDistributedBetaModel(basemodel,
                ncats)
        self.assertTrue(scipy.allclose(scipy.array([m.beta for m in
                gammamodel._models]), phydmslib.models.DiscreteGamma(
                gammamodel.alpha_lambda, gammamodel.beta_lambda,
                gammamodel.ncats)))
        for (param, pvalue) in paramvalues.items():
            if param != gammamodel.distributedparam:
                self.assertTrue(scipy.allclose(getattr(gammamodel, param),
                        pvalue))

        # try some updates and make sure everything remains OK
        for i in range(3):
            newvalues = {}
            for param in gammamodel.freeparams:
                (low, high) = gammamodel.PARAMLIMITS[param]
                if gammamodel.PARAMTYPES[param] == float:
                    newvalues[param] = random.uniform(low, high)
                else:
                    paramlength = gammamodel.PARAMTYPES[param][1]
                    newvalues[param] = scipy.random.uniform(
                            low, high, paramlength)
            gammamodel.updateParams(newvalues)
            self.assertTrue(scipy.allclose(scipy.array([m.beta for m in
                    gammamodel._models]), phydmslib.models.DiscreteGamma(
                    gammamodel.alpha_lambda, gammamodel.beta_lambda,
                    gammamodel.ncats)))
            for (param, pvalue) in newvalues.items():
                if param != gammamodel.distributedparam:
                    self.assertTrue(scipy.allclose(pvalue,
                            getattr(gammamodel, param)))
                    if param not in gammamodel.distributionparams:
                        self.assertTrue(all([scipy.allclose(pvalue,
                                getattr(m, param)) for m in
                                gammamodel._models]))

            # This is the opposite test of gammaomega
            self.assertTrue(gammamodel._models[0].branchScale >
                    gammamodel.branchScale >
                    gammamodel._models[-1].branchScale)

            t = 0.15
            for k in range(gammamodel.ncats):
                M = gammamodel.M(k, t)
                self.assertTrue(scipy.allclose(gammamodel._models[k].M(t), M))
                for param in gammamodel.freeparams:
                    if param not in gammamodel.distributionparams:
                        dM = gammamodel.dM(k, t, param, M)
                        self.assertTrue(scipy.allclose(dM,
                                gammamodel._models[k].dM(t, param, Mt=None)))

            # Check derivatives with respect to distribution params
            d_distparams = gammamodel.d_distributionparams
            self.assertTrue((d_distparams['alpha_lambda'] > 0).all())
            self.assertTrue((d_distparams['beta_lambda'] < 0).all())
            for param in gammamodel.distributionparams:
                diffs = []
                for k in range(gammamodel.ncats):
                    pvalue = getattr(gammamodel, param)
                    def func(x):
                        gammamodel.updateParams({param:x[0]})
                        return getattr(gammamodel._models[k],
                                gammamodel.distributedparam)
                    def dfunc(x):
                        gammamodel.updateParams({param:x[0]})
                        return gammamodel.d_distributionparams[param][k]
                    diff = scipy.optimize.check_grad(func, dfunc,
                            scipy.array([pvalue]))
                    gammamodel.updateParams({param:pvalue})
                    diffs.append(diff)
                diffs = scipy.array(diffs)
                self.assertTrue((diffs < 1e-5).all(), ("Excessive diff "
                        "for d_distributionparams[{0}] when "
                        "distributionparams = {1}:\n{2}".format(
                        param, gammamodel.distributionparams, diffs)))

            # Check the stationary state and deriviative of the stationary state
            # Stationary states should be different for each `k`
            self.assertFalse(all([scipy.allclose(gammamodel.stationarystate(i),
                    gammamodel.stationarystate(j)) for i in
                    range(gammamodel.ncats) for j in range(i+1, gammamodel.ncats)]))
            # The derviative of the stationary states should be different with
            # respect to beta for each `k`
            self.assertFalse(all([scipy.allclose(gammamodel.dstationarystate(i, "beta"),
                    gammamodel.dstationarystate(j, "beta")) for i in
                    range(gammamodel.ncats) for j in range(i+1, gammamodel.ncats)]))



class test_GammaDistributedBeta_ExpCM_empirical_phi(test_GammaDistributedBeta_ExpCM):
    """Test gamma distributed beta for `ExpCM_empirical_phi`."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    BASEMODEL = phydmslib.models.ExpCM_empirical_phi

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
