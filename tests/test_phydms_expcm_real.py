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
        """Tests command-line ``phydms`` with `ExpCM` on full data."""
#         prefs = '/Users/sarah/Google Drive/GS/lab/H1_H3_phyloComparison/H3_phydms/inputs/preferences/avgPrefsAll_H1H3noGaps.txt'
#         alignment = '/Users/sarah/Google Drive/GS/lab/H1_H3_phyloComparison/H3_phydms/inputs/sequences/H3_human_0_noStartnoStop_H1H3noGaps.fa'
#         tree = '/Users/sarah/Google Drive/GS/lab/H1_H3_phyloComparison/H3_phydms/inputs/sequences/RAxML_bestTree.H3_human_0'
        
        prefs = '/Users/sarah/Google Drive/GS/lab/H1_H3_phyloComparison/H3_phydms/inputs/preferences/HA_prefs_Doud2016.txt'
        alignment = '/Users/sarah/Google Drive/GS/lab/phydms_diversificationPressure/manuscript/phydms_diversifyingPressure/inputs/sequences/new/H1_HA_human_final_1.fasta'
        tree = '/Users/sarah/Google Drive/GS/lab/phydms_diversificationPressure/manuscript/phydms_diversifyingPressure/inputs/sequences/new/RAxML_bestTree.H1_HA_human_final_1'
        for f in [prefs, alignment]:
            self.assertTrue(os.path.isfile(f), "Can't find file {0}".format(f))
        outprefix = './phydmsReal_test_results/'
        if os.path.isdir(outprefix):
            shutil.rmtree(outprefix)

        subprocess.check_call(['phydms', alignment, tree, 
                'ExpCM_{0}'.format(prefs), outprefix])


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
