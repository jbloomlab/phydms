"""Tests `phydmslib.models.ExpCM_empirical_phi_divpressure` class.

Written by Jesse Bloom.
Edited by Sarah Hilton
"""


import random
import unittest
import scipy.linalg
import scipy.optimize
import numpy
from phydmslib.constants import (N_CODON, N_NT, AA_TO_INDEX, N_AA,
                                 CODON_NT_COUNT)
import phydmslib.models


class test_compare_ExpCM_emp_phi_with_without_divpress(unittest.TestCase):
    """Tests ``ExpCM`` with and without div pressure."""

    def test_compare(self):
        """Make sure all attributes are the same when `divpressure` is 0."""
        random.seed(1)
        numpy.random.seed(1)

        nsites = 6
        prefs = []
        minpref = 0.01
        for _r in range(nsites):
            rprefs = numpy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
        g = numpy.random.dirichlet([3] * N_NT)
        omega = 0.7
        omega2 = 0.2
        kappa = 2.5
        beta = 1.2
        divpressure = numpy.zeros(nsites)

        expcm = phydmslib.models.ExpCM_empirical_phi(
            prefs, g, omega=omega, kappa=kappa, beta=beta
        )

        expcm_divpressure = phydmslib.models.ExpCM_empirical_phi_divpressure(
            prefs,
            g,
            divPressureValues=divpressure,
            omega=omega,
            kappa=kappa,
            beta=beta,
            omega2=omega2,
        )

        self.assertTrue(numpy.allclose(expcm.stationarystate,
                                       expcm_divpressure.stationarystate),
                        "stationarystate differs.")
        self.assertTrue(numpy.allclose(expcm.Qxy, expcm_divpressure.Qxy),
                        "Qxy differs")
        self.assertTrue(
            numpy.allclose(expcm.Frxy, expcm_divpressure.Frxy), "Frxy differs")
        self.assertTrue(
            numpy.allclose(expcm.Prxy, expcm_divpressure.Prxy), "Prxy differs")
        t = 0.02
        self.assertTrue(
            numpy.allclose(expcm.M(t), expcm_divpressure.M(t)),
            "M({0}) differs".format(t))
        for param in ["kappa", "omega", "beta"]:
            self.assertTrue(
                numpy.allclose(
                    getattr(expcm, param), getattr(expcm_divpressure, param)),
                "param values differ for {0}".format(param))
            self.assertTrue(
                numpy.allclose(
                    expcm.dstationarystate(param),
                    (expcm_divpressure.dstationarystate(param))),
                "dstationarystate differs for {0}".format(param))
            self.assertTrue(
                numpy.allclose(
                    expcm.dM(t, param, expcm.M(t)),
                    (expcm_divpressure.dM(t, param, expcm_divpressure.M(t)))),
                "dM({0}) differs for {1}".format(t, param))


class testExpCM_empirical_phi_divpressure(unittest.TestCase):
    """Tests ``ExpCM`` with div pressure."""

    def test_ExpCM_empirical_phi_divpressure(self):
        """Init `ExpCM_empirical_phi_divpressure`, test, update, test again."""
        # create preferences
        random.seed(1)
        numpy.random.seed(1)
        self.nsites = 6
        self.prefs = []
        minpref = 0.01
        for _r in range(self.nsites):
            rprefs = numpy.random.dirichlet([0.5] * N_AA)
            rprefs[rprefs < minpref] = minpref
            rprefs /= rprefs.sum()
            self.prefs.append(dict(zip(sorted(AA_TO_INDEX.keys()), rprefs)))
        self.divpressure = numpy.random.randint(2, size=self.nsites)

        # create initial ExpCM
        g = numpy.random.dirichlet([3] * N_NT)
        omega = 0.7
        omega2 = 0.2
        kappa = 2.5
        beta = 1.2
        self.model = phydmslib.models.ExpCM_empirical_phi_divpressure(
            self.prefs,
            g=g,
            divPressureValues=self.divpressure,
            omega=omega,
            kappa=kappa,
            beta=beta,
            omega2=omega2)
        # now check ExpCM attributes / derivates, updating several times
        for _update in range(2):
            self.params = {
                "omega": random.uniform(0.1, 2),
                "kappa": random.uniform(0.5, 10),
                "beta": random.uniform(0.5, 3),
                "mu": random.uniform(0.05, 5.0),
                "omega2": random.uniform(0.1, 0.3),
            }
            self.model.updateParams(self.params)
            self.assertTrue(numpy.allclose(g, self.model.g))

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
                    nt_freqs[w] += self.model.prx[r][x] * CODON_NT_COUNT[w][x]
        self.assertTrue(numpy.allclose(sum(nt_freqs), 3 * self.nsites))
        nt_freqs = numpy.array(nt_freqs)
        nt_freqs /= nt_freqs.sum()
        self.assertTrue(
            numpy.allclose(nt_freqs, self.model.g, atol=1e-5),
            "Actual nt_freqs: {0}\nExpected (g): {1}"
            .format(nt_freqs, self.model.g))

        def func_phi(beta, expcm, w):
            expcm.updateParams({"beta": beta[0]})
            return expcm.phi[w]

        def func_dphi(beta, expcm, w):
            expcm.updateParams({"beta": beta[0]})
            return expcm.dphi_dbeta[w]

        for w in range(N_NT):
            diff = scipy.optimize.check_grad(func_phi,
                                             func_dphi,
                                             [self.model.beta],
                                             self.model,
                                             w,
                                             epsilon=1e-4)
            self.assertTrue(diff < 1e-4,
                            "dphi_dbeta diff {0} for w = {1}".format(diff, w))
        # back to original value
        self.model.updateParams(self.params)

    def check_dprx_dbeta(self):
        """Checks derivatives of `prx` with respect to `beta`."""

        def func_prx(beta, expcm, r, x):
            expcm.updateParams({"beta": beta[0]})
            return expcm.prx[r][x]

        def func_dprx(beta, expcm, r, x):
            expcm.updateParams({"beta": beta[0]})
            return expcm.dprx["beta"][r][x]

        for r in range(self.nsites):
            for x in range(N_CODON):
                diff = scipy.optimize.check_grad(func_prx,
                                                 func_dprx,
                                                 [self.model.beta],
                                                 self.model,
                                                 r,
                                                 x,
                                                 epsilon=1e-4,)
                self.assertTrue(diff < 1e-4,
                                "dprx_dbeta diff {0} for r = {1}, x = {2}"
                                .format(diff, r, x))
        # back to original value
        self.model.updateParams(self.params)

    def check_dQxy_dbeta(self):
        """Checks derivatives of `Qxy` with respect to `beta`."""

        def func_Qxy(beta, expcm, x, y):
            expcm.updateParams({"beta": beta[0]})
            return expcm.Qxy[x][y]

        def func_dQxy(beta, expcm, x, y):
            expcm.updateParams({"beta": beta[0]})
            return expcm.dQxy_dbeta[x][y]

        for x in random.sample(range(N_CODON), 3):
            for y in range(N_CODON):
                diff = scipy.optimize.check_grad(
                    func_Qxy,
                    func_dQxy,
                    [self.model.beta],
                    self.model,
                    x,
                    y,
                    epsilon=1e-4,
                )
                self.assertTrue(
                    diff < 1e-4,
                    "dQxy_dbeta diff {0} for x = {1}, y = {2}"
                    .format(diff, x, y))
        # back to original value
        self.model.updateParams(self.params)

    def check_ExpCM_attributes(self):
        """Make sure has the expected attribute values."""
        self.assertEqual(self.nsites, self.model.nsites)

        # make sure Prxy has rows summing to zero
        self.assertFalse(numpy.isnan(self.model.Prxy).any())
        self.assertFalse(numpy.isinf(self.model.Prxy).any())
        diag = numpy.eye(N_CODON, dtype="bool")
        for r in range(self.nsites):
            self.assertTrue(numpy.allclose(0, numpy.sum(self.model.Prxy[r],
                                           axis=1)))
            self.assertTrue(numpy.allclose(0, self.model.Prxy[r].sum()))
            self.assertTrue((self.model.Prxy[r][diag] <= 0).all())
            self.assertTrue((self.model.Prxy[r][~diag] >= 0).all())

        # make sure prx sums to 1 for each r
        self.assertTrue((self.model.prx >= 0).all())
        for r in range(self.nsites):
            self.assertTrue(numpy.allclose(1, self.model.prx[r].sum()))

        # prx is eigenvector or Prxy for the same r, but not different r
        for r in range(self.nsites):
            self.assertTrue(numpy.allclose(0, numpy.dot(self.model.prx[r],
                                                        self.model.Prxy[r])))
            if r > 0:
                (self.assertFalse(
                    numpy.allclose(0, numpy.dot(self.model.prx[r],
                                                self.model.Prxy[r - 1]))))

        # phi sums to one
        self.assertTrue(numpy.allclose(1, self.model.phi.sum()))

    def check_ExpCM_derivatives(self):
        """Makes sure derivatives are as expected."""
        # check derivatives of Prxy calculated by dPrxy
        # implementation looks a bit complex because `check_grad` function
        # can only be used for single values at a time, so have to loop
        # over r, x, y and so hash values to make faster

        def funcPrxy(paramvalue, paramname, expcm, r, x, y):
            if len(paramvalue) == 1:
                expcm.updateParams({paramname: paramvalue[0]})
            else:
                expcm.updateParams({paramname: paramvalue})
            return expcm.Prxy[r][x][y]

        def funcdPrxy(paramvalue, paramname, expcm, r, x, y):
            if len(paramvalue) == 1:
                expcm.updateParams({paramname: paramvalue[0]})
            else:
                expcm.updateParams({paramname: paramvalue})
            return expcm.dPrxy[paramname][r][x][y]

        for (pname, pvalue) in sorted(self.params.items())[::-1]:
            if pname == "mu":
                continue
            if isinstance(pvalue, float):
                pvalue = [pvalue]
            for r in random.sample(range(self.nsites), 2):  # check a few sites
                for x in random.sample(range(N_CODON), 3):  # check a few codon
                    for y in range(N_CODON):
                        diff = scipy.optimize.check_grad(
                            funcPrxy,
                            funcdPrxy,
                            pvalue,
                            pname,
                            self.model,
                            r,
                            x,
                            y,
                            epsilon=1e-4,)
                        (self.assertTrue(
                            diff < 1e-3,
                            "diff {0} for {1}: r = {2}, x = {3}, y = {4}, "
                            "beta = {5} pirAx = {6}, pirAy = {7}, mu = {8}, "
                            "omega = {9}, Frxy = {10}, Prxy = {11}, "
                            "phi = {12}, kappa = {13}, dQxy_dbeta = {14}, "
                            "dphi_dbeta = {15}, dPrxy_dbeta = {16}, "
                            "piAx_piAy_beta[r][x][y] = {17}".format(
                                    diff,
                                    pname,
                                    r,
                                    x,
                                    y,
                                    self.params["beta"],
                                    self.model.pi_codon[r][x],
                                    self.model.pi_codon[r][y],
                                    self.model.mu,
                                    self.model.omega,
                                    self.model.Frxy[r][x][y],
                                    self.model.Prxy[r][x][y],
                                    self.model.phi,
                                    self.model.kappa,
                                    self.model.dQxy_dbeta[x][y],
                                    self.model.dphi_dbeta,
                                    self.model.dPrxy["beta"][r][x][y],
                                    self.model.piAx_piAy_beta[r][x][y])))
            # back to original value
            self.model.updateParams(self.params)

    def check_ExpCM_matrix_exponentials(self):
        """Makes sure matrix exponentials are as expected."""
        for r in range(self.nsites):
            # fromdiag is recomputed Prxy after diagonalization
            fromdiag = numpy.dot(
                self.model.A[r],
                numpy.dot(numpy.diag(self.model.D[r]), self.model.Ainv[r]))
            self.assertTrue(
                numpy.allclose(self.model.Prxy[r], fromdiag, atol=1e-5),
                "Max diff {0}".format((self.model.Prxy[r] - fromdiag).max()))

            for t in [0.02, 0.2, 0.5]:
                direct = scipy.linalg.expm(self.model.Prxy[r] *
                                           self.model.mu * t)
                self.assertTrue(
                    numpy.allclose(self.model.M(t)[r], direct, atol=1e-6),
                    "Max diff {0}".format((self.model.M(t)[r] - direct).max()))
        # check derivatives of M calculated by dM
        # implementation looks a bit complex because `check_grad` function
        # can only be used for single values at a time, so have to loop
        # over r, x, y and so hash values to make faster

        def funcM(paramvalue, paramname, t, expcm, r, x, y, storedvalues):
            """Gets `M(t)[r][x][y]`."""
            key = ("M", tuple(paramvalue), paramname, t)
            if key not in storedvalues:
                if len(paramvalue) == 1:
                    expcm.updateParams({paramname: paramvalue[0]})
                else:
                    expcm.updateParams({paramname: paramvalue})
                storedvalues[key] = expcm.M(t)
            if len(storedvalues[key].shape) == 3:
                return storedvalues[key][r][x][y]
            else:
                return storedvalues[key][:, r, x, y]

        def funcdM(paramvalue, paramname, t, expcm, r, x, y, storedvalues):
            """Gets `dM`."""
            key = ("dM", tuple(paramvalue), paramname, t)
            if key not in storedvalues:
                if len(paramvalue) == 1:
                    expcm.updateParams({paramname: paramvalue[0]})
                else:
                    expcm.updateParams({paramname: paramvalue})
                storedvalues[key] = expcm.dM(t, paramname, expcm.M(t))
            if len(storedvalues[key].shape) == 3:
                return storedvalues[key][r][x][y]
            else:
                return storedvalues[key][:, r, x, y]

        for (pname, pvalue) in sorted(self.params.items()):
            storedvalues = {}  # used to hash values
            if isinstance(pvalue, float):
                pvalue = [pvalue]
            for t in [0.01, 0.2]:
                for r in random.sample(range(self.model.nsites), 2):
                    for x in random.sample(range(N_CODON), 3):
                        for y in range(N_CODON):
                            if x == y:
                                continue
                            diff = scipy.optimize.check_grad(
                                funcM,
                                funcdM,
                                pvalue,
                                pname,
                                t,
                                self.model,
                                r,
                                x,
                                y,
                                storedvalues,
                                epsilon=1e-4,
                            )
                            (
                                self.assertTrue(
                                    diff < 1e-3,
                                    "diff {0} for {1}: "
                                    "computed derivative = {10}, "
                                    "r = {2}, x = {3}, y = {4}, "
                                    "beta = {5}, pirAx = {6}, "
                                    "pirAy = {7}, t = {8}, mu = {9}".format(
                                        diff,
                                        pname,
                                        r,
                                        x,
                                        y,
                                        self.params["beta"],
                                        self.model.pi_codon[r][x],
                                        self.model.pi_codon[r][y],
                                        t,
                                        self.model.mu,
                                        funcdM(pvalue, pname, t, self.model,
                                               r, x, y, {}))))
                # back to original value
                self.model.updateParams(self.params)


if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
