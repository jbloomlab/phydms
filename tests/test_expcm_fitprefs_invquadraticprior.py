"""Tests invquadratic prior in `ExpCM_fitprefs`.

Written by Jesse Bloom."""


import random
import unittest
import copy
import scipy
import scipy.linalg
import sympy
from phydmslib.constants import *
import phydmslib.models


class test_ExpCM_fitprefs_invquadraticprior(unittest.TestCase):
    """Test `ExpCM_fitprefs` with inverse quadratic prior."""

    MODEL = phydmslib.models.ExpCM_fitprefs

    def setUp(self):
        """Set up for tests."""
        scipy.random.seed(1)
        random.seed(1)
        nsites = 1
        minpref = 0.001
        self.prefs = []
        for r in range(nsites):
            rprefs = scipy.random.dirichlet([0.7] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs[0] = rprefs[1] + 1.0e-8 # ensure near equal prefs handled OK
            rprefs /= rprefs.sum()
            self.prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
        self.expcm_fitprefs = self.MODEL(self.prefs, 
                prior=('invquadratic', 150.0, 0.5), kappa=3.0, omega=0.3,
                phi=scipy.random.dirichlet([5] * N_NT))
        assert len(self.expcm_fitprefs.zeta.flatten()) == nsites * (N_AA - 1)
        assert self.expcm_fitprefs.nsites == nsites


    def test_deriv(self):
        """Test analytical expression for derivative."""
        c1, c2, pi, theta = sympy.symbols('c1 c2 pi theta')
        lnpr = sympy.ln((1 / (1 + c1 * (pi - theta)**2))**c2)
        deriv = -2 * c1 * c2 * (pi - theta) / (1 + c1 * (pi - theta)**2)
        values = [(c1, 150), (c2, 0.5), (pi, 0.3), (theta, 0.421)]
        self.assertTrue(scipy.allclose(
                float(deriv.subs(values)), float(sympy.diff(lnpr, pi).subs(values))))


    def test_dlogprior(self):
        """Test `dlogprior`."""
        scipy.random.seed(1)

        expcm_fitprefs = copy.deepcopy(self.expcm_fitprefs)
        self.assertTrue(scipy.allclose(expcm_fitprefs.pi, expcm_fitprefs.origpi))
        if self.MODEL == phydmslib.models.ExpCM_fitprefs:
            newzeta = expcm_fitprefs.zeta.copy() * scipy.random.uniform(0.9, 1.0, 
                    expcm_fitprefs.zeta.shape)
        elif self.MODEL == phydmslib.models.ExpCM_fitprefs2:
            newzeta = expcm_fitprefs.zeta.copy() * scipy.random.uniform(0.01, 10.0, 
                    expcm_fitprefs.zeta.shape)
        else:
            raise RuntimeError("invalid MODEL: {0}".format(self.MODEL))
        expcm_fitprefs.updateParams({'zeta':newzeta})
        self.assertFalse(scipy.allclose(expcm_fitprefs.pi, expcm_fitprefs.origpi))
        nsites = expcm_fitprefs.nsites

        def func(zetari, i, r):
            zeta = expcm_fitprefs.zeta.copy()
            zeta.reshape(nsites, N_AA - 1)[r][i] = zetari
            expcm_fitprefs.updateParams({'zeta':zeta})
            return expcm_fitprefs.logprior

        def dfunc(zetari, i, r):
            zeta = expcm_fitprefs.zeta.copy()
            zeta.reshape(nsites, N_AA - 1)[r][i] = zetari
            expcm_fitprefs.updateParams({'zeta':zeta})
            return expcm_fitprefs.dlogprior('zeta')[i + r * (N_AA - 1)]

        j = 0
        for r in range(nsites):
            for i in range(N_AA - 1):
                zetari = scipy.array([expcm_fitprefs.zeta.reshape(
                        nsites, N_AA - 1)[r][i]])
                diff = scipy.optimize.check_grad(func, dfunc, zetari, i, r)
                deriv = dfunc(zetari, i, r)
                self.assertTrue(diff < max(1e-5, 1e-5 * abs(deriv)),
                        "{0}, {1}, {2}, {3}, {4}".format(
                        diff, zetari, i, r, deriv))
                j += 1


class test_ExpCM_fitprefs2_invquadraticprior(test_ExpCM_fitprefs_invquadraticprior):
    """Test `ExpCM_fitprefs2` with inverse quadratic prior."""

    MODEL = phydmslib.models.ExpCM_fitprefs2


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
