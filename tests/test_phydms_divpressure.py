"""Runs ``phydms`` using `ExpCM`.

This test examines the entire functionality of ``phydms``
when run from the command line using an `ExpCM` model with `--divpressure`.

Written by Jesse Bloom.
Edited by Sarah Hilton
"""

import os
import shutil
import unittest
import subprocess


class test_phydms_divpressure(unittest.TestCase):
    """Tests command-line ``phydms`` with `ExpCM` and `--divpressure`."""

    def test_NP_divpressure(self):
        """Tests command-line ``phydms`` with `ExpCM`  and `--divpressure` on simulated NP data."""
        dpdir = './divpressure_tests/'
        prefs = '{0}/NP_prefs.txt'.format(dpdir)
        alignment = '{0}/simulated_NP.fasta'.format(dpdir)
        tree = '{0}/NP_tree.newick'.format(dpdir)
        divpressure = '{0}/divpressure.txt'.format(dpdir)
        for f in [prefs, alignment]:
            self.assertTrue(os.path.isfile(f), "Can't find file {0}".format(f))
        outprefix = './NP_test_results/'
        if os.path.isdir(outprefix):
            shutil.rmtree(outprefix)

        subprocess.check_call(['phydms', alignment, tree, 
                'ExpCM_{0}'.format(prefs), outprefix, '--divpressure', divpressure])


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
