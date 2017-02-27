"""Tests `ExpCM` with fitting of preferences as free paramters.

Written by Jesse Bloom."""


import random
import unittest
import scipy
import scipy.linalg
import sympy
from phydmslib.constants import *
import phydmslib.models


class test_ExpCM_fitprefs(unittest.TestCase):
    """Test `ExpCM` with preferences as free parameters."""

    def test_DerivativeExpressions(self):
        """Makes sure we have right equations for derivatives."""
        random.seed(1)
        scipy.random.seed(1)
        omega, beta, pirAx, pirAy = sympy.symbols('omega beta pirAx pirAy')
        Frxy = omega * -beta * sympy.ln(pirAx / pirAy) / (1 - 
                (pirAx / pirAy)**beta)
        dFrxy_dpirAx = (-omega * beta / pirAx) * ((pirAx / pirAy)**beta * (
                sympy.ln((pirAx / pirAy)**beta) - 1) + 1) / ((1 - 
                (pirAx / pirAy)**beta)**2)
        dFrxy_dpirAx_prefsequal = -omega * beta / (2 * pirAx)
        dFrxy_dpirAy = (omega * beta / pirAy) * ((pirAx / pirAy)**beta * (
                sympy.ln((pirAx / pirAy)**beta) - 1) + 1) / ((1 - 
                (pirAx / pirAy)**beta)**2)
        dFrxy_dpirAy_prefsequal = omega * beta / (2 * pirAy)
        diffpref = 1.0e-5
        for itest in range(5):
            values = [[beta, random.uniform(0.5, 2.0)], 
                      [pirAx, random.uniform(0.01, 0.5)], 
                      [pirAy, random.uniform(0.01, 0.5)], 
                      [omega, random.uniform(0.1, 2.0)]]
            self.assertTrue(abs(values[1][1] - values[2][1]) > diffpref, 
                    "choose another random number seed as pirAx and pirAy "
                    "are too close.")
            self.assertTrue(scipy.allclose(float(dFrxy_dpirAx.subs(values)),
                    float(sympy.diff(Frxy, pirAx).subs(values))))
            self.assertTrue(scipy.allclose(float(dFrxy_dpirAy.subs(values)),
                    float(sympy.diff(Frxy, pirAy).subs(values))))
            values[2][1] = values[1][1] * (1 + diffpref)
            self.assertTrue(scipy.allclose(float(dFrxy_dpirAx_prefsequal.subs(
                    values)), float(sympy.diff(Frxy, pirAx).subs(values))))
            self.assertTrue(scipy.allclose(float(dFrxy_dpirAy_prefsequal.subs(
                    values)), float(sympy.diff(Frxy, pirAy).subs(values))))

    def test_ExpCM_fitprefs_derivs(self):
        """Initialize `ExpCM_fitprefs`, test derivatives with respect to `zeta`."""
        random.seed(1)
        scipy.random.seed(1)

        # initialize
        nsites = 1
        minpref = 0.001
        prefs = []
        for r in range(nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
        expcm_fitprefs = phydmslib.models.ExpCM_fitprefs(prefs, kappa=3.0,
                omega=0.3, mu=1.0, phi=scipy.random.dirichlet([5] * N_NT))
        assert len(expcm_fitprefs.zeta.flatten()) == nsites * (N_AA - 1)

        for r in range(nsites):
            for i in range(N_AA - 1):
                oldzeta = expcm_fitprefs.zeta.copy()
                oldpi = expcm_fitprefs.pi.copy()
                zeta = oldzeta.copy()
                zeta[r * nsites + i] *= 0.9
                expcm_fitprefs.updateParams({'zeta':zeta})
                self.assertFalse(scipy.allclose(oldzeta, expcm_fitprefs.zeta))
                self.assertFalse(scipy.allclose(oldpi[r], 
                        expcm_fitprefs.pi[r]))
                self.assertTrue(expcm_fitprefs.pi[r][i] > oldpi[r][i])
                self.assertTrue(all([expcm_fitprefs.pi[r][j] < oldpi[r][j]
                        for j in range(i + 1, N_AA)]))





if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
