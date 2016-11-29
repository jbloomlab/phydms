"""Tests `phydmslib.models.ExpCM_empirical_phi` class.

Written by Jesse Bloom.
"""


import random
import unittest
import scipy
import scipy.linalg
import sympy
from phydmslib.constants import *
import phydmslib.models


class testExpCM_empirical_phi(unittest.TestCase):

    def test_ExpCM_empirical_phi(self):
        """Initialize `ExpCM_empirical_phi`, test values, update, test again."""

        # create preferences
        random.seed(1)
        scipy.random.seed(1)
        self.nsites = 10
        self.prefs = []
        minpref = 0.01
        for r in range(self.nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref 
            rprefs /= rprefs.sum()
            self.prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))

        # create initial ExpCM
        g = scipy.random.dirichlet([3] * N_NT)
        omega = 0.7
        kappa = 2.5
        beta = 1.2
        self.expcm = phydmslib.models.ExpCM_empirical_phi(self.prefs, 
                g=g, omega=omega, kappa=kappa, beta=beta)
        self.assertTrue(scipy.allclose(g, self.expcm.g))

        # now check ExpCM attributes / derivates, updating several times
        for update in range(2):
            self.params = {'omega':random.uniform(0.1, 2),
                      'kappa':random.uniform(0.5, 10),
                      'beta':random.uniform(0.5, 5),
                      'mu':random.uniform(0.05, 5.0),
                     }
            self.expcm.updateParams(self.params)
            self.assertTrue(scipy.allclose(g, self.expcm.g))
            self.check_empirical_phi()
            self.check_dQxy_dbeta()
            self.check_dprx_dbeta()
            self.check_ExpCM_attributes()
            self.check_ExpCM_derivatives()
            self.check_ExpCM_matrix_exponentials()

    def check_empirical_phi(self):
        """Check that `phi` gives right `g`, and has right derivative."""
        nt_freqs = [0] * N_NT
        for r in range(self.nsites):
            for x in range(N_CODON):
                for w in range(N_NT):
                    nt_freqs[w] += self.expcm.prx[r][x] * CODON_NT_COUNT[w][x]
        self.assertTrue(scipy.allclose(sum(nt_freqs), 3 * self.nsites))
        nt_freqs = scipy.array(nt_freqs)
        nt_freqs /= nt_freqs.sum()
        self.assertTrue(scipy.allclose(nt_freqs, self.expcm.g, atol=1e-5), 
                "Actual nt_freqs: {0}\nExpected (g): {1}".format(
                nt_freqs, self.expcm.g))

        def func_phi(beta, expcm, w):
            expcm.updateParams({'beta':beta[0]})
            return expcm.phi[w]

        def func_dphi(beta, expcm, w):
            expcm.updateParams({'beta':beta[0]})
            return expcm.dphi_dbeta[w]

        for w in range(N_NT):
            diff = scipy.optimize.check_grad(func_phi, func_dphi, 
                    [self.expcm.beta], self.expcm, w, epsilon=1e-4)
            self.assertTrue(diff < 1e-4, 
                    "dphi_dbeta diff {0} for w = {1}".format(diff, w))
        self.expcm.updateParams(self.params) # back to original value

    def check_dprx_dbeta(self):
        """Checks derivatives of `prx` with respect to `beta`."""

        def func_prx(beta, expcm, r, x):
            expcm.updateParams({'beta':beta[0]})
            return expcm.prx[r][x]

        def func_dprx(beta, expcm, r, x):
            expcm.updateParams({'beta':beta[0]})
            return expcm.dprx['beta'][r][x]

        for r in range(self.nsites):
            for x in range(N_CODON):
                diff = scipy.optimize.check_grad(func_prx, func_dprx, 
                        [self.expcm.beta], self.expcm, r, x, epsilon=1e-4)
                self.assertTrue(diff < 1e-4, 
                        "dprx_dbeta diff {0} for r = {1}, x = {2}".format(
                        diff, r, x))
        self.expcm.updateParams(self.params) # back to original value

    def check_dQxy_dbeta(self):
        """Checks derivatives of `Qxy` with respect to `beta`."""

        def func_Qxy(beta, expcm, x, y):
            expcm.updateParams({'beta':beta[0]})
            return expcm.Qxy[x][y]

        def func_dQxy(beta, expcm, x, y):
            expcm.updateParams({'beta':beta[0]})
            return expcm.dQxy_dbeta[x][y]

        for x in random.sample(range(N_CODON), 3):
            for y in range(N_CODON):
                diff = scipy.optimize.check_grad(func_Qxy, func_dQxy, 
                        [self.expcm.beta], self.expcm, x, y, epsilon=1e-4)
                self.assertTrue(diff < 1e-4, 
                        "dQxy_dbeta diff {0} for x = {1}, y = {2}".format(
                        diff, x, y))
        self.expcm.updateParams(self.params) # back to original value

    def check_ExpCM_attributes(self):
        """Make sure has the expected attribute values."""
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

        # phi sums to one
        self.assertTrue(scipy.allclose(1, self.expcm.phi.sum()))

    def check_ExpCM_derivatives(self):
        """Makes sure derivatives are as expected."""
        # check derivatives of Prxy calculated by dPrxy
        # implementation looks a bit complex because `check_grad` function
        # can only be used for single values at a time, so have to loop 
        # over r, x, y and so hash values to make faster

        def funcPrxy(paramvalue, paramname, expcm, r, x, y):
            if len(paramvalue) == 1:
                expcm.updateParams({paramname:paramvalue[0]})
            else:
                expcm.updateParams({paramname:paramvalue})
            return expcm.Prxy[r][x][y]

        def funcdPrxy(paramvalue, paramname, expcm, r, x, y):
            if len(paramvalue) == 1:
                expcm.updateParams({paramname:paramvalue[0]})
            else:
                expcm.updateParams({paramname:paramvalue})
            return expcm.dPrxy[paramname][r][x][y]

        for (pname, pvalue) in sorted(self.params.items()):
            if pname == 'mu':
                continue
            if isinstance(pvalue, float):
                pvalue = [pvalue]
            for r in random.sample(range(self.nsites), 2): # check a few sites
                for x in random.sample(range(N_CODON), 3): # check a few codons
                    for y in range(N_CODON):
                        diff = scipy.optimize.check_grad(funcPrxy, funcdPrxy, 
                                pvalue, pname, self.expcm, r, x, y, epsilon=1e-4)
                        self.assertTrue(diff < 1e-3, ("diff {0} for {1}:" +
                                " r = {2}, x = {3}, y = {4}, beta = {5} " +
                                "pirAx = {6}, pirAy = {7}, mu = {8}, " +
                                "omega = {9}, Frxy = {10}, Prxy = {11}, " +
                                "phi = {12}, kappa = {13}, dQxy_dbeta = {14}, " +
                                "dphi_dbeta = {15}, dPrxy_dbeta = {16}"
                                ).format(diff, pname, r, x, y, 
                                self.params['beta'], self.expcm.pi_codon[r][x], 
                                self.expcm.pi_codon[r][y], self.expcm.mu,
                                self.expcm.omega, self.expcm.Frxy[r][x][y],
                                self.expcm.Prxy[r][x][y], self.expcm.phi, 
                                self.expcm.kappa, self.expcm.dQxy_dbeta[x][y],
                                self.expcm.dphi_dbeta, 
                                self.expcm.dPrxy['beta'][r][x][y]))
            self.expcm.updateParams(self.params) # back to original value

    def check_ExpCM_matrix_exponentials(self):
        """Makes sure matrix exponentials are as expected."""
        for r in range(self.nsites):
            # fromdiag is recomputed Prxy after diagonalization
            fromdiag = scipy.dot(self.expcm.A[r], scipy.dot(scipy.diag(
                    self.expcm.D[r]), self.expcm.Ainv[r]))
            self.assertTrue(scipy.allclose(self.expcm.Prxy[r], fromdiag,
                    atol=1e-5), "Max diff {0}".format(
                    (self.expcm.Prxy[r] - fromdiag).max()))

            for t in [0.02, 0.2, 0.5]:
                direct = scipy.linalg.expm(self.expcm.Prxy[r] * self.expcm.mu * t)
                self.assertTrue(scipy.allclose(self.expcm.M(t)[r], direct, atol=1e-6),
                        "Max diff {0}".format((self.expcm.M(t)[r] - direct).max()))
        # check derivatives of M calculated by dM
        # implementation looks a bit complex because `check_grad` function
        # can only be used for single values at a time, so have to loop 
        # over r, x, y and so hash values to make faster
        def funcM(paramvalue, paramname, t, expcm, r, x, y, storedvalues):
            """Gets `M(t)[r][x][y]`."""
            key = ('M', tuple(paramvalue), paramname, t)
            if key not in storedvalues:
                if len(paramvalue) == 1:
                    expcm.updateParams({paramname:paramvalue[0]})
                else:
                    expcm.updateParams({paramname:paramvalue})
                storedvalues[key] = expcm.M(t)
            if len(storedvalues[key].shape) == 3:
                return storedvalues[key][r][x][y]
            else:
                return storedvalues[key][:,r,x,y]
        def funcdM(paramvalue, paramname, t, expcm, r, x, y, storedvalues):
            """Gets `dM`."""
            key = ('dM', tuple(paramvalue), paramname, t)
            if key not in storedvalues:
                if len(paramvalue) == 1:
                    expcm.updateParams({paramname:paramvalue[0]})
                else:
                    expcm.updateParams({paramname:paramvalue})
                storedvalues[key] = expcm.dM(t, paramname, expcm.M(t))
            if len(storedvalues[key].shape) == 3:
                return storedvalues[key][r][x][y]
            else:
                return storedvalues[key][:,r,x,y]
        for (pname, pvalue) in sorted(self.params.items()):
            storedvalues = {} # used to hash values
            if isinstance(pvalue, float):
                pvalue = [pvalue]
            for t in [0.01, 0.2]:
                for r in random.sample(range(self.expcm.nsites), 2):
                    for x in random.sample(range(N_CODON), 3):
                        for y in range(N_CODON):
                            if x == y:
                                continue
                            diff = scipy.optimize.check_grad(funcM, funcdM, pvalue,
                                    pname, t, self.expcm, r, x, y, storedvalues,
                                    epsilon=1e-4) 
                            self.assertTrue(diff < 1e-3, ("diff {0} for {1}:" +
                                " computed derivative = {10}, " +
                                " r = {2}, x = {3}, y = {4}, beta = {5}, " +
                                "pirAx = {6}, pirAy = {7}, t = {8}, mu = {9}"
                                ).format(diff, pname, r, x, y, 
                                self.params['beta'], self.expcm.pi_codon[r][x], 
                                self.expcm.pi_codon[r][y], t, self.expcm.mu,
                                funcdM(pvalue, pname, t, self.expcm, 
                                r, x, y, {})))
                self.expcm.updateParams(self.params) # back to original value



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
