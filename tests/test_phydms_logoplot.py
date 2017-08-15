"""Runs ``phydms_logoplot``

This test makes sure ``phydms_logoplot`` runs OK. It is not
able to test the actual output for correctness is this is a plot.

Written by Jesse Bloom.
"""

import os
import unittest
import subprocess


class test_phydms_logoplot(unittest.TestCase):
    """Tests command-line ``phydms_logoplot``."""

    def setUp(self):
        """Define input data."""
        self.diffprefs = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './expected_NP_test_results/ExpCM_NP_prefs_diffprefsbysite.txt'))
        self.omega = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './expected_NP_test_results/ExpCM_NP_prefs_omegabysite.txt'))
        self.prefs = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './NP_data/NP_prefs.tsv'))
        self.outprefix = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './logoplot_test_results/'))
        self.stringency = '2.99'
        if not os.path.isdir(self.outprefix):
            os.mkdir(self.outprefix)

    def test_prefs_logoplot(self):
        """Tests ``--prefs`` option to ``phydms_logoplot``."""
        plotfile = os.path.join(self.outprefix,
                'prefs_logoplot.pdf')
        if os.path.isfile(plotfile):
            os.remove(plotfile)
        subprocess.call(['phydms_logoplot', '--prefs', self.prefs,
                plotfile, '--stringency', self.stringency,
                '--nperline', '72', '--mapmetric', 'charge'])
        self.assertTrue(os.path.isfile(plotfile))

    def test_diffprefs_logoplot(self):
        """Tests ``--diffprefs`` option to ``phydms_logoplot``."""
        plotfile = os.path.join(self.outprefix,
                'diffprefs_logoplot.pdf')
        if os.path.isfile(plotfile):
            os.remove(plotfile)
        subprocess.call(['phydms_logoplot', '--diffprefs', self.diffprefs,
                plotfile, '--nperline', '72'])
        print(plotfile)
        self.assertTrue(os.path.isfile(plotfile))

    def test_omegaoverlay(self):
        """Tests ``--omegabysite`` option to ``phydms_logoplot``."""
        plotfile = os.path.join(self.outprefix,
                'omegaoverlay_logoplot.pdf')
        if os.path.isfile(plotfile):
            os.remove(plotfile)
        subprocess.call(['phydms_logoplot', '--prefs', self.prefs,
                plotfile, '--stringency', self.stringency,
                '--nperline', '72', '--omegabysite', self.omega,
                '--minP', '0.001'])
        self.assertTrue(os.path.isfile(plotfile))



if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
