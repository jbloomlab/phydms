"""Runs the model adequacy protocol.

This test examines the functionality of the model adequcy protocol and checks
the results against the output in `expected_modeladequacy_results`.

The test tests the function with `ncpus` is equal to 1 (no multiprocessing)
and with `ncpus` greater than 1 (multiprocessing).

Written by Sarah Hilton.
"""

import unittest
import phydmslib.simulate
import phydmslib.file_io
import phydmslib.modeladequacy
from phydmslib.constants import *
import os
from statsmodels.sandbox.stats.multicomp import multipletests
import scipy
import numpy as np
import itertools
import pandas as pd
import cProfile
import pstats
import glob
import subprocess


class test_modeladequacy_ExpCM_mp(unittest.TestCase):
    """Runs model adequacy on an ExpCM with >1 CPU."""
    # run parameters
    MODEL = "ExpCM_{0}".format(os.path.join(os.path.abspath(os.path.dirname(__file__)),
                               "modeladequacy_tests/HA_short_prefs_nsites10.csv"))
    EXPECTED = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                            "expected_modeladequacy_results/ExpCM_pvalues_seed0_2000rep_pvalues.csv")
    SEED = 0
    TREE = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        "modeladequacy_tests/HA_short_nsites10_nseqs34_tree.newick")
    NCPUS = 4

    def test_modeladequacy(self):
        """Runs model adequacy and compares against expected results."""
        """Runs model adequacy and compares against expected results."""
        n_sim = 2000
        alignment = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                 "modeladequacy_tests/HA_short_nsites10_nseqs34.fasta")
        outprefix = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                 "_model_adequacy_results")
        cmd = ["phydms_modeladequacy", outprefix, alignment,
               self.MODEL, "--number_simulations", str(n_sim), "--tree",
               self.TREE, "--seed", str(self.SEED), "--ncpus", str(self.NCPUS)]
        subprocess.check_call(cmd)

        final = (pd.read_csv("{0}_pvalues.csv".format(outprefix))
                 .sort_values(by=["site", "metric"]))
        expected = (pd.read_csv(self.EXPECTED)
                    .sort_values(by=["site", "metric"]))

        self.assertTrue(scipy.allclose(final["pvalue"], expected["pvalue"]),
                       " pvalue: Expected \n{0}\n \nvs.\n \n{1}.".format(
                       final["pvalue"].to_string(),
                       final["pvalue"].to_string()))
        self.assertTrue(scipy.allclose(final["qvalue"], expected["qvalue"]),
                       " qvalue: Expected \n{0}\n \nvs.\n \n{1}.".format(
                       final["qvalue"].to_string(),
                       final["qvalue"].to_string()))

        # remove files
        for fname in glob.glob("{0}_*".format(outprefix)):
            os.remove(fname)


class test_modeladequacy_ExpCM_noMP(test_modeladequacy_ExpCM_mp):
    """Runs model adequacy on an ExpCM with 1 CPU."""
    # run parameters
    NCPUS = 1


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
