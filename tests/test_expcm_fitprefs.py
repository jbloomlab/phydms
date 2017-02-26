"""Tests `ExpCM` with fitting of preferences as free paramters.

Written by Jesse Bloom."""


import random
import unittest
import scipy
import scipy.linalg
import sympy
from phydmslib.constants import *
import phydmslib.models


class testExpCM_fitprefs(unittest.TestCase):
    """Test `ExpCM` with preferences as free parameters."""

    def test_DerivativeExpressions(self):
        """Makes sure we have right equations for derivatives."""
        random.seed(1)
        omega, beta, pirAx, pirAy = sympy.symbols('omega beta pirAx pirAy')
        Frxy = omega * -beta * sympy.ln(pirAx / pirAy) / (1 - 
                (pirAx / pirAy)**beta)
        dFrxy_dpirAx = omega * beta * ((pirAx / pirAy)**beta * (1 - sympy.ln(
                (pirAx / pirAy)**beta)) - 1) / (pirAx * (1 - 
                (pirAx / pirAy)**beta)**2)
        dFrxy_dpirAx_prefsequal = -omega * beta / (2 * pirAx)
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
            values[2][1] = values[1][1] * (1 + diffpref)
            self.assertTrue(scipy.allclose(float(dFrxy_dpirAx_prefsequal.subs(
                    values)), float(sympy.diff(Frxy, pirAx).subs(values))))



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
