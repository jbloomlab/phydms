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
        tree = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './divpressure_tests/NP_tree.newick'))
        alignment = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './divpressure_tests/simulated_NP.fasta'))
        prefs = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './divpressure_tests/NP_prefs.txt'))
        divpressure = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './divpressure_tests/divpressure.txt'))
        divpressure2 = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './divpressure_tests/divpressure2.csv'))
        for f in [prefs, alignment]:
            self.assertTrue(os.path.isfile(f), "Can't find file {0}".format(f))
        outprefix = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './divpressure_test_results/'))
        if outprefix[-1] != "/":
            outprefix = "{0}/".format(outprefix)
        if os.path.isdir(outprefix):
            shutil.rmtree(outprefix)

        
        subprocess.check_call(['phydms_divpressure', outprefix, alignment,
                prefs, "--tree", tree, divpressure, divpressure2])
        actual = os.path.abspath(os.path.join(outprefix,
                "modelcomparison.csv"))
        actual = pd.read_csv(actual)
        expected = os.path.abspath(os.path.join(os.path.dirname(__file__),
                "./expected_divpressure_test_results/modelcomparison.csv"))
        expected = pd.read_csv(expected)
        assert((actual.columns.values == expected.columns.values).all())
        columns = [x for x in actual.columns.values if x != "value"]
        actual.sort_values(by=columns, inplace=True)
        expected.sort_values(by=columns, inplace=True)
        self.assertTrue(scipy.allclose(actual["value"],
                expected["value"], atol=1e-2, rtol=1e-5))



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
