"""Tests `phydmslib.models.ExpCM` class.

Written by Jesse Bloom.

Uses `sympy` to make sure attributes and derivatives of attributes
are correct for `ExpCM` implemented in `phydmslib.models`."""


import random
import unittest
import scipy
import scipy.linalg
import sympy
from phydmslib.constants import *
import phydmslib.models


class testExpCM(unittest.TestCase):

    def test_ExpCM(self):
        """Initialize `ExpCM`, test values, update, test again."""

        # create preferences
        random.seed(1)
        scipy.random.seed(1)
        self.nsites = 2
        self.prefs = []
        minpref = 0.01
        for r in range(self.nsites):
            rprefs = scipy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref 
            rprefs /= rprefs.sum()
            self.prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))

        # create initial ExpCM
        phi = scipy.random.dirichlet([2] * N_NT)
        omega = 0.7
        kappa = 2.5
        beta = 1.9
        self.expcm = phydmslib.models.ExpCM(self.prefs, phi=phi, omega=omega,
                kappa=kappa, beta=beta)
        self.assertTrue(scipy.allclose(phi, self.expcm.phi))
        self.assertTrue(scipy.allclose(omega, self.expcm.omega))
        self.assertTrue(scipy.allclose(kappa, self.expcm.kappa))
        self.assertTrue(scipy.allclose(beta, self.expcm.beta))

        # now check ExpCM attributes / derivates, updating several times
        for update in range(2):
            self.params = {'omega':random.uniform(*self.expcm.PARAMLIMITS['omega']),
                      'kappa':random.uniform(*self.expcm.PARAMLIMITS['kappa']),
                      'beta':random.uniform(0.5, 2.5),
                      'eta':scipy.array([random.uniform(*self.expcm.PARAMLIMITS['eta']) for i in range(N_NT - 1)]),
                      'mu':random.uniform(0.05, 3.0),
                     }
            self.expcm.updateParams(self.params)
            self.check_ExpCM_attributes()
            self.check_ExpCM_derivatives()
            self.check_ExpCM_matrix_exponentials()

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

    def check_ExpCM_derivatives(self):
        """Use `sympy` to check values and derivatives of `ExpCM` attributes."""
        (Prxy, Qxy, phiw, beta, omega, eta0, eta1, eta2, kappa) = sympy.symbols(
                'Prxy, Qxy, phiw, beta, omega, eta0, eta1, eta2, kappa')

        values = {'beta':self.params['beta'],
                  'omega':self.params['omega'],
                  'kappa':self.params['kappa'],
                  'eta0':self.params['eta'][0],
                  'eta1':self.params['eta'][1],
                  'eta2':self.params['eta'][2],
                  }

        # check Prxy
        for r in range(self.nsites):
            for x in range(N_CODON):
                pirAx = self.prefs[r][INDEX_TO_AA[CODON_TO_AA[x]]]
                for y in [yy for yy in range(N_CODON) if yy != x]:
                    pirAy = self.prefs[r][INDEX_TO_AA[CODON_TO_AA[y]]]
                    if not CODON_SINGLEMUT[x][y]:
                        Prxy = 0
                    else:
                        w = NT_TO_INDEX[[ynt for (xnt, ynt) in zip(INDEX_TO_CODON[x], 
                                INDEX_TO_CODON[y]) if xnt != ynt][0]]
                        if w == 0:
                            phiw = 1 - eta0                            
                        elif w == 1:
                            phiw = eta0 * (1 - eta1)
                        elif w == 2:
                            phiw = eta0 * eta1 * (1 - eta2)
                        elif w == 3:
                            phiw = eta0 * eta1 * eta2
                        else:
                            raise ValueError("Invalid w")
                        self.assertTrue(scipy.allclose(float(phiw.subs(values)), 
                                self.expcm.phi[w]))
                        if CODON_TRANSITION[x][y]:
                            Qxy = kappa * phiw
                        else:
                            Qxy = phiw
                        self.assertTrue(scipy.allclose(float(Qxy.subs(values)), 
                                self.expcm.Qxy[x][y]))
                        if CODON_NONSYN[x][y]:
                            if pirAx == pirAy:
                                Prxy = Qxy * omega
                            else:
                                Prxy = Qxy * omega * (-beta * scipy.log(pirAx / pirAy) / (1 - (pirAx / pirAy)**beta))
                        else:
                            Prxy = Qxy
                    for (name, actual, expect) in [
                            ('Prxy', self.expcm.Prxy[r][x][y], Prxy),
                            ('dPrxy_dkappa', self.expcm.dPrxy['kappa'][r][x][y], sympy.diff(Prxy, kappa)),
                            ('dPrxy_domega', self.expcm.dPrxy['omega'][r][x][y], sympy.diff(Prxy, omega)),
                            ('dPrxy_dbeta', self.expcm.dPrxy['beta'][r][x][y], sympy.diff(Prxy, beta)),
                            ('dPrxy_deta0', self.expcm.dPrxy['eta'][0][r][x][y], sympy.diff(Prxy, eta0)),
                            ('dPrxy_deta1', self.expcm.dPrxy['eta'][1][r][x][y], sympy.diff(Prxy, eta1)),
                            ('dPrxy_deta2', self.expcm.dPrxy['eta'][2][r][x][y], sympy.diff(Prxy, eta2)),
                            ]:
                        if Prxy == 0:
                            expectval = 0
                        else:
                            expectval = float(expect.subs(values))
                        self.assertTrue(scipy.allclose(actual, expectval, atol=1e-4),
                                "{0}: {1} vs {2}".format(name, actual, expectval))

        # check prx
        qxs = [sympy.Symbol('qx{0}'.format(x)) for x in range(N_CODON)]
        frxs = [sympy.Symbol('frx{0}'.format(x)) for x in range(N_CODON)]
        prx = sympy.Symbol('prx')
        phixs = [sympy.Symbol('phix{0}'.format(w)) for w in range(3)]
        for r in range(self.nsites):
            for x in range(N_CODON):
                pirAx = self.prefs[r][INDEX_TO_AA[CODON_TO_AA[x]]]
                frxs[x] = pirAx**beta
                xcodon = INDEX_TO_CODON[x]
                assert len(phixs) == len(xcodon)
                for (w, xwnt) in enumerate(xcodon):
                    xw = NT_TO_INDEX[xwnt]
                    if xw == 0:
                        phixs[w] = 1 - eta0
                    elif xw == 1:
                        phixs[w] = eta0 * (1 - eta1)
                    elif xw == 2:
                        phixs[w] = eta0 * eta1 * (1 - eta2)
                    elif xw == 3:
                        phixs[w] = eta0 * eta1 * eta2
                    else:
                        raise ValueError("invalid xw")
                qxs[x] = phixs[0] * phixs[1] * phixs[2]
            for x in range(N_CODON):
                prx = frxs[x] * qxs[x] / sum(frx * qx for (frx, qx) in zip(frxs, qxs))
                for (name, actual, expect) in [
                        ('prx', self.expcm.prx[r][x], prx),
                        ('dprx_dbeta', self.expcm.dprx['beta'][r][x], sympy.diff(prx, beta)),
                        ('dprx_deta0', self.expcm.dprx['eta'][0][r][x], sympy.diff(prx, eta0)),
                        ('dprx_deta1', self.expcm.dprx['eta'][1][r][x], sympy.diff(prx, eta1)),
                        ('dprx_deta2', self.expcm.dprx['eta'][2][r][x], sympy.diff(prx, eta2)),
                        ]:
                    expectval = float(expect.subs(values))
                    self.assertTrue(scipy.allclose(actual, expectval, atol=1e-5),
                            "{0}: {1} vs {2}".format(name, actual, expectval))

    def check_ExpCM_matrix_exponentials(self):
        """Makes sure matrix exponentials of ExpCM are as expected."""
        for r in range(self.nsites):
            # fromdiag is recomputed Prxy after diagonalization
            fromdiag = scipy.dot(self.expcm.A[r], scipy.dot(scipy.diag(
                    self.expcm.D[r]), self.expcm.Ainv[r]))
            self.assertTrue(scipy.allclose(self.expcm.Prxy[r], fromdiag,
                    atol=1e-5), "Max diff {0}".format((self.expcm.Prxy[r] - fromdiag).max()))

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
            for t in [0.01, 0.2, 0.5]:
                for r in range(self.expcm.nsites):
                    for x in range(N_CODON):
                        for y in range(N_CODON):
                            diff = scipy.optimize.check_grad(funcM, funcdM, pvalue,
                                    pname, t, self.expcm, r, x, y, storedvalues) 
                            self.assertTrue(diff < 2e-3, ("diff {0} for {1}:" +
                                " r = {2}, x = {3}, y = {4}, beta = {5} " +
                                "pirAx = {6}, pirAy = {7}, t = {8}, mu = {9}"
                                ).format(diff, pname, r, x, y, 
                                self.params['beta'], self.expcm.pi_codon[r][x], 
                                self.expcm.pi_codon[r][y], t, self.expcm.mu))
                self.expcm.updateParams(self.params) # back to original value



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
