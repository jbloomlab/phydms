"""Tests ``phydms`` flag --divpressure on a simulated dataset of NP sequences under diversifying selection.

Written by Sarah Hilton."""


import sys
import os
import unittest
import subprocess
import dms_tools.file_io



class TestOnDivPresureNPs(unittest.TestCase):
    """
    Runs ``phydms`` --divpressure on test data, compares to known results.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.expected_dir = './expected_divpressure_results/'
        self.test_dir = './divpressure_test_results/'
        self.prefs = 'shortNP_prefs.txt'
        self.alignment = 'simulatedNPs.fasta'
        self.divpressure = 'divpressure.txt'
        self.tree = "NP_tree"
        for f in [self.prefs, self.alignment, self.divpressure, self.tree]:
            self.assertTrue(os.path.isfile(f), "Can't find required file {0}".format(f))
        self.likelihood_files = 'divpressure_loglikelihood.txt'
        self.params_files = 'divpressure_modelparams.txt'
        for f in [self.likelihood_files, self.params_files]:
            fname = "%s/%s" % (self.expected_dir, f)
            self.assertTrue(os.path.isfile(fname), "Cannot find required file %s" % fname)
            toremove = "%s/%s" % (self.test_dir, f)
            if os.path.isfile(toremove):
                os.remove(toremove)
# 
# 
    def test_divpressure(self):
        """Runs ``phydms`` with --divpressure flag."""
        cmds = ['phydms', self.alignment, self.tree, 'ExpCM_'+self.prefs, self.test_dir+"divpressure", '--divpressure', self.divpressure]
        subprocess.call(cmds)
# 
        sys.stderr.write('\nTesting for presence of expected output files...\n')
        for f in [self.likelihood_files,self.params_files]:
            fname = '%s/%s' % (self.test_dir, f)
            self.assertTrue(os.path.isfile(fname), "Failed to created expected file %s" % fname)
# 
        sys.stderr.write('\nTesting for expected likelihood values...\n')
        for f in [self.likelihood_files]:
            with open('%s/%s' % (self.expected_dir, f)) as fin:
                expected = float(fin.read().split('=')[-1])
            with open('%s/%s' % (self.test_dir, f)) as fin:
                actual = float(fin.read().split('=')[-1])
            self.assertTrue(abs(expected - actual) < 1, "Unexpectedly large differences in %s: %g versus %g" % (f, expected, actual))
#
# 
        sys.stderr.write('\nTesting for expected parameter values...\n')
        for f in [self.params_files]:
            with open('%s/%s' % (self.expected_dir, f)) as fin:
                expected = dict([(line.split('=')[0], float(line.split('=')[1])) for line in fin.readlines()])
            with open('%s/%s' % (self.test_dir, f)) as fin:
                actual = dict([(line.split('=')[0], float(line.split('=')[1])) for line in fin.readlines()])
            self.assertTrue(set(expected.keys()) == set(actual.keys()), "Different parameters in %s: %s versus %s" % (f, str(expected.keys()), str(actual.keys())))
            for key in actual.keys():
                self.assertTrue(abs(actual[key] - expected[key]) < 0.02, "Unexpectedly large difference for %s in %s: %g versus %g" % (key, f, actual[key], expected[key]))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
