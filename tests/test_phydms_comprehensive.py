"""Runs ``phydms_comprehensive``

This test examines the entire functionality of ``phydms_comprehensive``
when run from the command line.

Written by Jesse Bloom and Sarah Hilton.
"""

import os
import shutil
import unittest
import subprocess


class test_phydms_comprehensive(unittest.TestCase):
    """Tests command-line ``phydms_comprehensive``."""

    def test_NP(self):
        """Tests command-line ``phydms_comprehensive`` on NP data."""
        npdir = './NP_data/'
        prefs = '{0}/NP_prefs.tsv'.format(npdir)
        alignment = '{0}/NP_alignment.fasta'.format(npdir)
        tree = '{0}/NP_tree.newick'.format(npdir)
        for f in [prefs, alignment]:
            self.assertTrue(os.path.isfile(f), "Can't find file {0}".format(f))
        outprefix = './NP_test_results/'

        subprocess.check_call(['phydms_comprehensive', outprefix, alignment,
                prefs, "--tree", tree, "--omegabysite", '--brlen', 'scale'])


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
