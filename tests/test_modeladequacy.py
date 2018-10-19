"""Tests `phydmslib.modeladequacy` class.

Written by Sarah Hilton.
"""


import random
import unittest
import scipy
from phydmslib.constants import *
import phydmslib.models
import phydmslib.modeladequacy
import pandas as pd


class test_modeladequacy(unittest.TestCase):

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

        self.assertTrue(scipy.allclose(scipy.repeat(1.0, self.nsites),
                                       self.expcm.stationarystate.sum(axis=1)))

        # other tests
        self.check_stationarystate_aminoacid_frequencies()

    def check_stationarystate_aminoacid_frequencies(self):
        """Make sure the function `make_stationary_state_prefs` is correct."""
        amino_acids = list(INDEX_TO_AA.values())
        amino_acids.sort()
        header = ["site"] + [INDEX_TO_AA[x] for x in range(20)]
        df = []
        for r in range(self.expcm.nsites):
            aa_ss = scipy.zeros(N_AA)
            for x in range(N_CODON):
                aa = CODON_TO_AA[x]
                aa_ss[aa] += self.expcm.stationarystate[r][x]
            assert scipy.isclose(aa_ss.sum(), 1.0)
            df.append([r+1]+list(aa_ss))
        df = pd.DataFrame(df, columns=header)

        self.ss_aa_freqs = (phydmslib.modeladequacy
                            .calc_stationary_state_freqs(self.expcm))

        # check that the stationary state frequencies are the same
        self.assertTrue(scipy.allclose(self.ss_aa_freqs,
                                       df[amino_acids].values))

        # check that the pref sets are the same
        self.ss_aa_prefs = (phydmslib.modeladequacy
                            .make_stationary_state_prefs(self.expcm))
        self.assertTrue(df.equals(self.ss_aa_prefs))
        self.assertTrue(self.ss_aa_prefs.equals(df))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
