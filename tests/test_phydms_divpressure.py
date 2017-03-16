"""Runs ``phydms_testdivpressure``

This test examines the entire functionality of ``phydms_testdivpressure``
when run from the command line.

Written by Jesse Bloom and Sarah Hilton.
"""

import os
import shutil
import unittest
import subprocess
import pandas as pd
import scipy


class test_phydms_testdivpressure(unittest.TestCase):
    """Tests command-line ``phydms_testdivpressure``."""

    def test_NP(self):
        """Tests command-line ``phydms_testdivpressure`` on NP data."""
        npdir = './divpressure_tests/'
        prefs = '{0}/NP_prefs.txt'.format(npdir)
        alignment = '{0}/simulated_NP.fasta'.format(npdir)
        divpressure = '{0}/divpressure.txt'.format(npdir)
        divpressure2 = '{0}/divpressure2.csv'.format(npdir)
        tree = '{0}/NP_tree.newick'.format(npdir)
        for f in [prefs, alignment]:
            self.assertTrue(os.path.isfile(f), "Can't find file {0}".format(f))
        outprefix = './divpressure_test_results/'
        if os.path.isdir(outprefix):
            shutil.rmtree(outprefix)


        subprocess.check_call(['phydms_divpressure', outprefix, alignment,
                prefs, "--tree", tree, divpressure, divpressure2])

        actual = pd.read_csv(outprefix + "modelcomparison.csv")
        expected = pd.read_csv("./expected_divpressure_test_results/" \
                + "modelcomparison.csv")
        assert((actual.columns.values == expected.columns.values).all())
        columns = [x for x in actual.columns.values if x != "value"]
        actual.sort_values(by=columns, inplace=True)
        expected.sort_values(by=columns, inplace=True)
        self.assertTrue(scipy.allclose(actual["value"],
                expected["value"], atol=1e-2, rtol=1e-5))



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
