"""Tests `phydmslib.models.GammaDistributedOmegaModel` class.

Written by Jesse Bloom.
"""


import random
import unittest
import numpy
import scipy
import scipy.linalg
import scipy.optimize
from phydmslib.constants import *
import phydmslib.models


class test_GammaDistributedOmega_ExpCM(unittest.TestCase):
    """Test gamma distributed omega for `ExpCM`."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    BASEMODEL = phydmslib.models.ExpCM

    def test_GammaDistributedOmega(self):
        """Initialize, test values, update, test again."""

        random.seed(1)
        scipy.random.seed(1)
        nsites = 10

        if self.BASEMODEL == phydmslib.models.ExpCM:
            prefs = []
            minpref = 0.01
            for r in range(nsites):
                rprefs = scipy.random.dirichlet([0.5] * N_AA)
                rprefs[rprefs < minpref] = minpref
                rprefs /= rprefs.sum()
                prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
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
        elif self.BASEMODEL == phydmslib.models.YNGKP_M0:
            e_pw = scipy.random.uniform(0.4, 0.6, size=(3, N_NT))
            e_pw = e_pw / e_pw.sum(axis=1, keepdims=True)
            basemodel = self.BASEMODEL(e_pw, nsites)
            paramvalues = {
                        'kappa':2.5,
                        'omega':0.7,
                        'mu':0.5,
                        }
            assert set(paramvalues.keys()) == set(basemodel.freeparams)
            basemodel.updateParams(paramvalues)
        else:
            raise ValueError("Invalid BASEMODEL: {0}".format(self.BASEMODEL))
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))

        ncats = 4
        gammamodel = phydmslib.models.GammaDistributedOmegaModel(basemodel,
                ncats)
        self.assertTrue(numpy.allclose(numpy.array([m.omega for m in
                gammamodel._models]), phydmslib.models.DiscreteGamma(
                gammamodel.alpha_lambda, gammamodel.beta_lambda,
                gammamodel.ncats)))
        for (param, pvalue) in paramvalues.items():
            if param != gammamodel.distributedparam:
                self.assertTrue(numpy.allclose(getattr(gammamodel, param),
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
            self.assertTrue(numpy.allclose(numpy.array([m.omega for m in
                    gammamodel._models]), phydmslib.models.DiscreteGamma(
                    gammamodel.alpha_lambda, gammamodel.beta_lambda,
                    gammamodel.ncats)))
            for (param, pvalue) in newvalues.items():
                if param != gammamodel.distributedparam:
                    self.assertTrue(numpy.allclose(pvalue,
                            getattr(gammamodel, param)))
                    if param not in gammamodel.distributionparams:
                        self.assertTrue(all([numpy.allclose(pvalue,
                                getattr(m, param)) for m in
                                gammamodel._models]))

            self.assertTrue(gammamodel._models[0].branchScale <
                    gammamodel.branchScale <
                    gammamodel._models[-1].branchScale)

            t = 0.15
            for k in range(gammamodel.ncats):
                M = gammamodel.M(k, t)
                self.assertTrue(numpy.allclose(gammamodel._models[k].M(t), M))
                for param in gammamodel.freeparams:
                    if param not in gammamodel.distributionparams:
                        dM = gammamodel.dM(k, t, param, M)
                        self.assertTrue(numpy.allclose(dM,
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
                            numpy.array([pvalue]))
                    gammamodel.updateParams({param:pvalue})
                    diffs.append(diff)
                diffs = numpy.array(diffs)
                self.assertTrue((diffs < 1e-5).all(), ("Excessive diff "
                        "for d_distributionparams[{0}] when "
                        "distributionparams = {1}:\n{2}".format(
                        param, gammamodel.distributionparams, diffs)))



class test_GammaDistributedOmega_YNGKP_M0(test_GammaDistributedOmega_ExpCM):
    """Test gamma distributed omega for `YNGKP_M0` (this is `YNGKP_M5`)."""

    # use approach here to run multiple tests:
    # http://stackoverflow.com/questions/17260469/instantiate-python-unittest-testcase-with-arguments
    BASEMODEL = phydmslib.models.YNGKP_M0

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
