"""Tests `phydmslib.models.YNGKP_M0` class.

Written by Jesse Bloom.
Edited by Sarah Hilton

Uses `sympy` to make sure attributes and derivatives of attributes
are correct for `YNGKP_M0` implemented in `phydmslib.models`."""


import random
import unittest
import scipy
import scipy.linalg
import sympy
from phydmslib.constants import *
import phydmslib.models


class testYNGKP_M0(unittest.TestCase):

    def test_YNGKP_M0(self):
        """Initialize `YNGKP_M0`, test values, update, test again."""

        # create alignment
        sequences = [""]
        #self.e_pa = scipy.array([[4.0, 0.0, 0.0, 0.0], [0.0, 1.0, 2.0, 1.0], [0.0, 2.0, 0.0, 2.0]])
        scipy.random.seed(1)
        random.seed(1)
        self.e_pa = scipy.random.uniform(0.12, 1.0, size = (3, N_NT))
        self.e_pa = self.e_pa/self.e_pa.sum(axis=1, keepdims=True)
        self.nsites = 10

        # create initial YNGKP_M0
        phi = scipy.random.dirichlet([2] * N_NT)
        omega = 0.4
        kappa = 2.5
        self.YNGKP_M0 = phydmslib.models.YNGKP_M0(self.e_pa, self.nsites, omega=omega, kappa=kappa)
        self.assertTrue(scipy.allclose(omega, self.YNGKP_M0.omega))
        self.assertTrue(scipy.allclose(kappa, self.YNGKP_M0.kappa))

        # now check YNGKP_M0 attributes / derivates, updating several times
        for update in range(2):
            self.params = {'omega':random.uniform(0.1, 2),
                      'kappa':random.uniform(0.5, 10),
                      'mu':random.uniform(0.05, 5.0),
                     }
            self.YNGKP_M0.updateParams(self.params)
            self.check_YNGKP_M0_attributes()
            self.check_YNGKP_M0_derivatives()
            self.check_YNGKP_M0_matrix_exponentials()

    def check_YNGKP_M0_attributes(self):
        """Make sure `YNGKP_M0` has the expected attribute values."""
        self.assertEqual(self.nsites, self.YNGKP_M0.nsites)

        # make sure Prxy has rows summing to zero
        self.assertFalse(scipy.isnan(self.YNGKP_M0.Pxy).any())
        self.assertFalse(scipy.isinf(self.YNGKP_M0.Pxy).any())
        diag = scipy.eye(N_CODON, dtype='bool')
        self.assertTrue(scipy.allclose(0, scipy.sum(self.YNGKP_M0.Pxy[0],
                axis=1)))
        self.assertTrue(scipy.allclose(0, self.YNGKP_M0.Pxy[0].sum()))
        self.assertTrue((self.YNGKP_M0.Pxy[0][diag] <= 0).all())
        self.assertTrue((self.YNGKP_M0.Pxy[0][~diag] >= 0).all())
        assert self.YNGKP_M0.Pxy.shape == (1, N_CODON, N_CODON)
        for param in self.YNGKP_M0.dPxy:
            assert self.YNGKP_M0.dPxy[param].shape == (1, N_CODON, N_CODON)
            assert self.YNGKP_M0.B[param].shape == (1, N_CODON, N_CODON)
        assert self.YNGKP_M0.Ainv.shape == (1, N_CODON, N_CODON)
        assert self.YNGKP_M0.A.shape == (1, N_CODON, N_CODON)
        assert self.YNGKP_M0.M(0.2).shape == (self.nsites, N_CODON, N_CODON)
        for (pname, value) in self.params.items():
            assert self.YNGKP_M0.dM(0.2, pname, self.YNGKP_M0.M(0.2)).shape == (self.nsites, N_CODON, N_CODON)


    def check_YNGKP_M0_derivatives(self):
        """Makes sure derivatives are as expected."""
        # check derivatives of Prxy calculated by dPrxy
        # implementation looks a bit complex because `check_grad` function
        # can only be used for single values at a time, so have to loop
        # over r, x, y and so hash values to make faster

        def funcPxy(paramvalue, paramname, YNGKP_M0, x, y):
            if len(paramvalue) == 1:
                YNGKP_M0.updateParams({paramname:paramvalue[0]})
            else:
                YNGKP_M0.updateParams({paramname:paramvalue})
            return YNGKP_M0.Pxy[0][x][y]

        def funcdPxy(paramvalue, paramname, YNGKP_M0, x, y):
            if len(paramvalue) == 1:
                YNGKP_M0.updateParams({paramname:paramvalue[0]})
            else:
                YNGKP_M0.updateParams({paramname:paramvalue})
            return YNGKP_M0.dPxy[paramname][0][x][y]

        for (pname, pvalue) in sorted(self.params.items())[::-1]:
            if pname == 'mu':
                continue
            if isinstance(pvalue, float):
                pvalue = [pvalue]
            for x in random.sample(range(N_CODON), 3): # check a few codons
                for y in range(N_CODON):
                    diff = scipy.optimize.check_grad(funcPxy, funcdPxy,
                            pvalue, pname, self.YNGKP_M0, x, y, epsilon=1e-4)
                    self.assertTrue(diff < 1e-3, ("diff {0} for {1}:" +
                            "x = {2}, y = {3} " +
                            "mu = {4}, " +
                            "omega = {5}, Pxy = {6}, " +
                            "kappa = {7}"
                            ).format(diff, pname, x, y,
                            self.YNGKP_M0.mu, self.YNGKP_M0.omega,
                            self.YNGKP_M0.Pxy[0][x][y],
                            self.YNGKP_M0.kappa))
            self.YNGKP_M0.updateParams(self.params) # back to original value

    def check_YNGKP_M0_matrix_exponentials(self):
        """Makes sure matrix exponentials of YNGKP_M0 are as expected."""
        for r in range(1):
            # fromdiag is recomputed Prxy after diagonalization
            fromdiag = scipy.dot(self.YNGKP_M0.A[0], scipy.dot(scipy.diag(
                    self.YNGKP_M0.D[r]), self.YNGKP_M0.Ainv[r]))
            self.assertTrue(scipy.allclose(self.YNGKP_M0.Pxy[r], fromdiag,
                    atol=1e-5), "Max diff {0}".format((self.YNGKP_M0.Pxy[r] - fromdiag).max()))

            for t in [0.02, 0.2, 0.5]:
                direct = scipy.linalg.expm(self.YNGKP_M0.Pxy[r] * self.YNGKP_M0.mu * t)
                self.assertTrue(scipy.allclose(self.YNGKP_M0.M(t)[r], direct, atol=1e-6),
                        "Max diff {0}".format((self.YNGKP_M0.M(t)[r] - direct).max()))
        # check derivatives of M calculated by dM
        # implementation looks a bit complex because `check_grad` function
        # can only be used for single values at a time, so have to loop
        # over r, x, y and so hash values to make faster
        def funcM(paramvalue, paramname, t, YNGKP_M0, r, x, y, storedvalues):
            """Gets `M(t)[r][x][y]`."""
            key = ('M', tuple(paramvalue), paramname, t)
            if key not in storedvalues:
                if len(paramvalue) == 1:
                    YNGKP_M0.updateParams({paramname:paramvalue[0]})
                else:
                    YNGKP_M0.updateParams({paramname:paramvalue})
                storedvalues[key] = YNGKP_M0.M(t)
            if len(storedvalues[key].shape) == 3:
                return storedvalues[key][r][x][y]
            else:
                return storedvalues[key][:,r,x,y]
        def funcdM(paramvalue, paramname, t, YNGKP_M0, r, x, y, storedvalues):
            """Gets `dM`."""
            key = ('dM', tuple(paramvalue), paramname, t)
            if key not in storedvalues:
                if len(paramvalue) == 1:
                    YNGKP_M0.updateParams({paramname:paramvalue[0]})
                else:
                    YNGKP_M0.updateParams({paramname:paramvalue})
                storedvalues[key] = YNGKP_M0.dM(t, paramname, YNGKP_M0.M(t))
            if len(storedvalues[key].shape) == 3:
                return storedvalues[key][r][x][y]
            else:
                return storedvalues[key][:,r,x,y]
        for (pname, pvalue) in sorted(self.params.items()):
            storedvalues = {} # used to hash values
            if isinstance(pvalue, float):
                pvalue = [pvalue]
            for t in [0.01, 0.2, 0.5]:
                for r in range(1):
                    for x in range(N_CODON):
                        for y in range(N_CODON):
                            diff = scipy.optimize.check_grad(funcM, funcdM, pvalue,
                                    pname, t, self.YNGKP_M0, r, x, y, storedvalues)
                            self.assertTrue(diff < 1e-3, ("diff {0} for {1}:" +
                                " r = {2}, x = {3}, y = {4}"
                                ).format(diff, pname, r, x, y))
                self.YNGKP_M0.updateParams(self.params) # back to original value



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
