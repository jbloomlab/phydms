"""Runs ``phydms`` using `YNGKP_M0`.

This test examines the entire functionality of ``phydms``
when run from the command line using an `YNGKP_M0` model.

Written by Jesse Bloom.
Edited by Sarah Hilton
"""

import os
import shutil
import unittest
import subprocess


class test_phydms_divpressure(unittest.TestCase):
    """Tests command-line ``phydms`` with `YNGKP_M0`."""

    def test_NP_YNGKP_M0(self):
        """Tests command-line ``phydms`` with `ExpCM`  and `--divpressure` on simulated NP data."""
        dpdir = './NP_data/'
        alignment = '{0}/NP_alignment.fasta'.format(dpdir)
        #alignment = '{0}/short_fasta_test.fa'.format(dpdir)
        tree = '{0}/NP_tree.newick'.format(dpdir)
        self.assertTrue(os.path.isfile(alignment), "Can't find file {0}".format(alignment))
        outprefix = './NP_test_results/'
        if os.path.isdir(outprefix):
            shutil.rmtree(outprefix)

        subprocess.check_call(['phydms', alignment, tree,
                'YNGKP_M0', outprefix,])


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
