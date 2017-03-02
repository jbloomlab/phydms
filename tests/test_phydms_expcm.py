"""Runs ``phydms`` using `ExpCM`.

This test examines the entire functionality of ``phydms``
when run from the command line using an `ExpCM` model.

Written by Jesse Bloom.
"""

import os
import shutil
import unittest
import subprocess


class test_phydms_ExpCM(unittest.TestCase):
    """Tests command-line ``phydms`` with `ExpCM`."""

    def test_NP(self):
        """Tests command-line ``phydms`` with `ExpCM` on NP data."""
        npdir = './NP_data/'
        prefs = '{0}/NP_prefs.tsv'.format(npdir)
        alignment = '{0}/NP_alignment.fasta'.format(npdir)
        tree = '{0}/NP_tree.newick'.format(npdir)
        for f in [prefs, alignment]:
            self.assertTrue(os.path.isfile(f), "Can't find file {0}".format(f))
        outprefix = './NP_test_results/ExpCM'
        if os.path.isdir(outprefix):
            shutil.rmtree(outprefix)

        subprocess.check_call(['phydms', alignment, tree, 
                'ExpCM_{0}'.format(prefs), outprefix, '--omegabysite', '--diffprefsbysite'])


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
