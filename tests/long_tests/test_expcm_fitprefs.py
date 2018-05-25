"""Tests `ExpCM` with fitting of preferences as free parameters.

Written by Jesse Bloom."""


import random
import unittest
import copy
import scipy
import scipy.linalg
import sympy
from phydmslib.constants import *
import phydmslib.models


class test_ExpCM_fitprefs(unittest.TestCase):
    """Test `ExpCM_fitprefs`."""

    MODEL = phydmslib.models.ExpCM_fitprefs

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
            values = [[beta, 1], 
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

        expcm_fitprefs = copy.deepcopy(self.expcm_fitprefs)
        for r in range(expcm_fitprefs.nsites):
            for x in range(N_CODON):
                for y in range(N_CODON):
                    if x == y:
                        continue
                    values = {}
                    values[beta] = expcm_fitprefs.beta
                    values[omega] = expcm_fitprefs.omega
                    values[pirAx] = expcm_fitprefs.pi_codon[r][x]
                    values[pirAy] = expcm_fitprefs.pi_codon[r][y]
                    Qxy = expcm_fitprefs.Qxy[x][y]

                    # check Prxy values
                    if values[pirAx] == values[pirAy]:
                        if CODON_TO_AA[x] == CODON_TO_AA[y]:
                            self.assertTrue(scipy.allclose(Qxy,
                                    expcm_fitprefs.Prxy[r][x][y]))
                        else:
                            self.assertTrue(scipy.allclose(Qxy * values[omega],
                                    expcm_fitprefs.Prxy[r][x][y]))
                    else:
                        self.assertTrue(scipy.allclose(Qxy * float(Frxy.subs(
                            values.items())), expcm_fitprefs.Prxy[r][x][y]))

                    # check dFrxy_dpi
                    if values[pirAx] == values[pirAy]:
                        if CODON_TO_AA[x] == CODON_TO_AA[y]:
                            self.assertTrue(scipy.allclose(0,
                                    -expcm_fitprefs.tildeFrxy[r][x][y] / 
                                    values[pirAx]))
                            self.assertTrue(scipy.allclose(0,
                                    -expcm_fitprefs.tildeFrxy[r][x][y] / 
                                    values[pirAy]))
                        else:
                            self.assertTrue(scipy.allclose(
                                    float(dFrxy_dpirAx_prefsequal.subs(
                                    values.items())),
                                    -expcm_fitprefs.tildeFrxy[r][x][y] / 
                                    values[pirAx]))
                            self.assertTrue(scipy.allclose(
                                    float(dFrxy_dpirAy_prefsequal.subs(
                                    values.items())),
                                    expcm_fitprefs.tildeFrxy[r][x][y] / 
                                    values[pirAy]))
                    else:
                        self.assertTrue(scipy.allclose(
                                float(dFrxy_dpirAx.subs(values.items())),
                                -expcm_fitprefs.tildeFrxy[r][x][y] / 
                                values[pirAx]))
                        self.assertTrue(scipy.allclose(
                                float(dFrxy_dpirAy.subs(values.items())),
                                expcm_fitprefs.tildeFrxy[r][x][y] / 
                                values[pirAy]))

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
                prior=None, kappa=3.0, omega=0.3,
                phi=scipy.random.dirichlet([5] * N_NT))
        assert len(self.expcm_fitprefs.zeta.flatten()) == nsites * (N_AA - 1)
        assert self.expcm_fitprefs.nsites == nsites


    def test_origbeta(self):
        """Test `origbeta` parameter."""
        for origbeta in [0.8, 1.0, 1.2]:
            expcm_fitprefs2 = self.MODEL(self.prefs,
                    prior=None,
                    kappa=self.expcm_fitprefs.kappa,
                    omega=self.expcm_fitprefs.omega,
                    mu=self.expcm_fitprefs.mu,
                    phi=self.expcm_fitprefs.phi,
                    origbeta=origbeta)
            zeta = self.expcm_fitprefs.zeta
            zeta2 = expcm_fitprefs2.zeta
            pi = self.expcm_fitprefs.zeta
            pi2 = expcm_fitprefs2.zeta
            self.assertTrue(scipy.allclose(zeta, zeta2) == (origbeta == 1))
            self.assertTrue(scipy.allclose(pi, pi2) == (origbeta == 1))
            for r in range(expcm_fitprefs2.nsites):
                for x in range(N_CODON):
                    pix = self.expcm_fitprefs.pi_codon[r][x]
                    pix2 = expcm_fitprefs2.pi_codon[r][x]
                    for y in range(N_CODON):
                        piy = self.expcm_fitprefs.pi_codon[r][y]
                        piy2 = expcm_fitprefs2.pi_codon[r][y]
                        if piy == pix:
                            self.assertTrue(pix2 == piy2)
                        elif origbeta == 1:
                            self.assertTrue((pix == pix2) and (piy == piy2))
                        elif piy > pix:
                            self.assertTrue((origbeta > 1) == (piy2 / pix2 >
                                    piy / pix))
                        else:
                            self.assertTrue((origbeta > 1) == (piy2 / pix2 <
                                    piy / pix))


    def test_zeta_updates(self):
        """Test updating `zeta`."""
        random.seed(1)
        scipy.random.seed(1)

        expcm_fitprefs = copy.deepcopy(self.expcm_fitprefs)

        self.assertTrue(scipy.allclose(expcm_fitprefs.pi, expcm_fitprefs.origpi))
        k = 0
        for r in range(expcm_fitprefs.nsites):
            for i in range(N_AA - 1):
                oldzeta = expcm_fitprefs.zeta.copy()
                oldpi = expcm_fitprefs.pi.copy()
                zeta = oldzeta.copy()
                zeta[k] *= 0.9
                expcm_fitprefs.updateParams({'zeta':zeta})
                self.assertFalse(scipy.allclose(oldzeta, expcm_fitprefs.zeta))
                self.assertFalse(scipy.allclose(oldpi[r], 
                        expcm_fitprefs.pi[r]))
                self.assertFalse(scipy.allclose(expcm_fitprefs.pi,
                        expcm_fitprefs.origpi))
                if self.MODEL == phydmslib.models.ExpCM_fitprefs:
                    self.assertTrue(expcm_fitprefs.pi[r][i] > oldpi[r][i])
                    self.assertTrue(all([expcm_fitprefs.pi[r][j] < 
                            oldpi[r][j] for j in range(i + 1, N_AA)]))
                elif self.MODEL == phydmslib.models.ExpCM_fitprefs2:
                    self.assertTrue(expcm_fitprefs.pi[r][i] < oldpi[r][i])
                    zeta[k] *= 2
                    expcm_fitprefs.updateParams({'zeta':zeta})
                    self.assertTrue(expcm_fitprefs.pi[r][i] > oldpi[r][i])
                k += 1

    def test_dprx_dzeta(self):
        """Test `dprx['zeta']`."""
        random.seed(1)

        expcm_fitprefs = copy.deepcopy(self.expcm_fitprefs)
        nsites = expcm_fitprefs.nsites

        def func(zetari, i, r, x):
            zeta = expcm_fitprefs.zeta.copy()
            zeta.reshape(nsites, N_AA - 1)[r][i] = zetari
            expcm_fitprefs.updateParams({'zeta':zeta})
            return expcm_fitprefs.prx[r][x]

        def dfunc(zetari, i, r, x):
            zeta = expcm_fitprefs.zeta.copy()
            zeta.reshape(nsites, N_AA - 1)[r][i] = zetari
            expcm_fitprefs.updateParams({'zeta':zeta})
            return expcm_fitprefs.dprx['zeta'][i + r * (N_AA - 1)][r][x]

        j = 0
        for r in range(nsites):
            for i in range(N_AA - 1):
                zetari = scipy.array([expcm_fitprefs.zeta.reshape(
                        nsites, N_AA - 1)[r][i]])
                for x in range(N_CODON):
                    diff = scipy.optimize.check_grad(func, dfunc,
                            zetari, i, r, x)
                    deriv = dfunc(zetari, i, r, x)
                    self.assertTrue(diff < max(1e-4, 1e-5 * abs(deriv)),
                            "{0}, {1}, {2}, {3}, {4}, {5}, {6}".format(
                            diff, zetari, i, r, x, CODON_TO_AA[x], deriv))
                j += 1


    def test_dPrxy_dzeta(self):
        """Test `dPrxy['zeta']`."""
        random.seed(1)

        expcm_fitprefs = copy.deepcopy(self.expcm_fitprefs)
        nsites = expcm_fitprefs.nsites

        def func(zetari, i, r, x, y):
            zeta = expcm_fitprefs.zeta.copy()
            zeta.reshape(nsites, N_AA - 1)[r][i] = zetari
            expcm_fitprefs.updateParams({'zeta':zeta})
            return expcm_fitprefs.Prxy[r][x][y]

        def dfunc(zetari, i, r, x, y):
            zeta = expcm_fitprefs.zeta.copy()
            zeta.reshape(nsites, N_AA - 1)[r][i] = zetari
            expcm_fitprefs.updateParams({'zeta':zeta})
            return expcm_fitprefs.dPrxy['zeta'][i + r * (N_AA - 1)][r][x][y]

        j = 0
        for r in range(nsites):
            for i in range(N_AA - 1):
                zetari = scipy.array([expcm_fitprefs.zeta.reshape(
                        nsites, N_AA - 1)[r][i]])
                for x in random.sample(range(N_CODON), 10):
                    for y in random.sample(range(N_CODON), 10):
                        if x == y:
                            continue
                        diff = scipy.optimize.check_grad(func, dfunc,
                                zetari, i, r, x, y)
                        deriv = expcm_fitprefs.dPrxy['zeta'][j][r][x][y]
                        self.assertTrue(diff < max(1e-4, 1e-5 * abs(deriv)),
                                "{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}".format(
                                diff, zetari, i, r, x, y, CODON_TO_AA[x], 
                                CODON_TO_AA[y], deriv))
                j += 1

    def test_dM_dzeta(self):
        """Test `dM['zeta']`."""
        random.seed(1)

        expcm_fitprefs = copy.deepcopy(self.expcm_fitprefs)
        nsites = expcm_fitprefs.nsites

        def func(zetari, i, r, x, y, t):
            zeta = expcm_fitprefs.zeta.copy()
            zeta.reshape(nsites, N_AA - 1)[r][i] = zetari
            expcm_fitprefs.updateParams({'zeta':zeta})
            return expcm_fitprefs.M(t)[r][x][y]

        def dfunc(zetari, i, r, x, y, t):
            zeta = expcm_fitprefs.zeta.copy()
            zeta.reshape(nsites, N_AA - 1)[r][i] = zetari
            expcm_fitprefs.updateParams({'zeta':zeta})
            j = i + r * (N_AA - 1)
            return expcm_fitprefs.dM(t, 'zeta', None)[j][r][x][y]

        for r in range(nsites):
            for i in range(N_AA - 1):
                zetari = scipy.array([expcm_fitprefs.zeta.reshape(
                        nsites, N_AA - 1)[r][i]])
                for x in random.sample(range(N_CODON), 10):
                    for y in random.sample(range(N_CODON), 10):
                        if x == y:
                            continue
                        for t in [0.1, 0.5]:
                            diff = scipy.optimize.check_grad(func, dfunc,
                                    zetari, i, r, x, y, t)
                            deriv = dfunc(zetari, i, r, x, y, t)
                            self.assertTrue(diff < max(0.02, 1e-4 * abs(deriv)),
                                    "{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, "
                                    "{8}, {9}".format(
                                    diff, zetari, i, r, x, y, CODON_TO_AA[x], 
                                    CODON_TO_AA[y], deriv, t))


class test_ExpCM_fitprefs2(test_ExpCM_fitprefs):
    """Test `ExpCM_fitprefs2`."""

    MODEL = phydmslib.models.ExpCM_fitprefs2

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
