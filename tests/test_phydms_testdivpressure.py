"""Runs ``phydms_testdivpressure``

This test examines the entire functionality of ``phydms_testdivpressure``
when run from the command line.

Written by Jesse Bloom and Sarah Hilton.
"""

import os
import shutil
import unittest
import subprocess


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
        outprefix = './divpressure_tests/results/'
        if os.path.isdir(outprefix):
            shutil.rmtree(outprefix)


        subprocess.check_call(['phydms_testdivpressure', outprefix, alignment,
                prefs, "--tree", tree, divpressure, divpressure2])


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
