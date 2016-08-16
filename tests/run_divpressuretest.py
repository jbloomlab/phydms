"""Tests ``phydms`` flag --divpressure on a simulated dataset of NP sequences under diversifying selection.

Written by Sarah Hilton."""


import sys
import os
import random
import unittest
import subprocess
import dms_tools.file_io



class TestOnDivPresureNPs(unittest.TestCase):
    """
    Runs ``phydms`` --divpressure on test data, compares to known results.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        random.seed(1)
        self.expected_dir = './expected_divpressure_results/'
        self.test_dir = './divpressure_test_results/'
        self.prefs = self.expected_dir + 'shortNP_prefs.txt'
        self.alignment = self.expected_dir + 'simulatedNPs.fasta'
        self.divpressure = self.expected_dir + 'divpressure.txt'
        with open(self.divpressure) as f:
            lines = [line.split() for line in f.readlines() if not line.isspace() and line[0] != '#']
        assert all([len(line) == 2 for line in lines])
        sites = [tup[0] for tup in lines]
        divpressures = [tup[1] for tup in lines]
        self.neg_divpressure = self.expected_dir + 'neg_divpressure.txt'
        with open(self.neg_divpressure, 'w') as f:
            for (site, dp) in zip(sites, divpressures):
                dp = str(-2.0 * float(dp))
                f.write('{0} {1}\n'.format(site, dp))
        random.shuffle(divpressures)
        self.rand_divpressure = self.expected_dir + 'rand_divpressure.txt'
        with open(self.rand_divpressure, 'w') as f:
            f.write('\n'.join('{0} {1}'.format(site, dp) for (site, dp) in zip(sites, divpressures)))
        self.tree = self.expected_dir + "NP_tree.newick"
        for f in [self.prefs, self.alignment, self.divpressure, self.rand_divpressure, self.tree]:
            self.assertTrue(os.path.isfile(f), "Can't find required file {0}".format(f))
        self.likelihood_suffix = 'loglikelihood.txt'
        self.params_suffix = 'modelparams.txt'
        self.prefixes = ['divpressure', 'no_divpressure', 'rand_divpressure', 'neg_divpressure']
        for suffix in [self.likelihood_suffix, self.params_suffix]:
            for prefix in self.prefixes:
                fname = "{0}/{1}_{2}".format(self.expected_dir, prefix, suffix)
                self.assertTrue(os.path.isfile(fname), "Cannot find required file %s" % fname)
                toremove = "{0}/{1}_{2}".format(self.test_dir, prefix, suffix)
            if os.path.isfile(toremove):
                os.remove(toremove)
 
 
    def test_divpressure(self):
        """Runs ``phydms`` with --divpressure flag."""
        for cmds in [
                ['phydms', self.alignment, self.tree, 'ExpCM_'+self.prefs, self.test_dir+"divpressure", '--divpressure', self.divpressure],
                ['phydms', self.alignment, self.tree, 'ExpCM_'+self.prefs, self.test_dir+"neg_divpressure", '--divpressure', self.neg_divpressure],
                ['phydms', self.alignment, self.tree, 'ExpCM_'+self.prefs, self.test_dir+"no_divpressure"],
                ['phydms', self.alignment, self.tree, 'ExpCM_'+self.prefs, self.test_dir+"rand_divpressure", '--divpressure', self.rand_divpressure]
                ]:
            subprocess.check_call(cmds)
 
        sys.stderr.write('\nTesting for presence of expected output files...\n')
 
        sys.stderr.write('\nTesting for expected likelihood and parameter values...\n')
        for prefix in self.prefixes:
            # test for likelihood
            f_expected = '{0}/{1}_{2}'.format(self.expected_dir, prefix, self.likelihood_suffix)
            f_actual = '{0}/{1}_{2}'.format(self.test_dir, prefix, self.likelihood_suffix)
            with open(f_expected) as fin:
                expected = float(fin.read().split('=')[-1])
            with open(f_actual) as fin:
                actual = float(fin.read().split('=')[-1])
            self.assertTrue(abs(expected - actual) < 1, "Unexpectedly large differences in %s: %g versus %g" % (f_actual, expected, actual))
            # test for expected parameter values
            f_expected = '{0}/{1}_{2}'.format(self.expected_dir, prefix, self.params_suffix)
            f_actual = '{0}/{1}_{2}'.format(self.test_dir, prefix, self.params_suffix)
            with open(f_expected) as fin:
                expected = dict([(line.split('=')[0], float(line.split('=')[1])) for line in fin.readlines()])
            with open(f_actual) as fin:
                actual = dict([(line.split('=')[0], float(line.split('=')[1])) for line in fin.readlines()])
            self.assertTrue(set(expected.keys()) == set(actual.keys()), "Different parameters in %s: %s versus %s" % (f_actual, str(expected.keys()), str(actual.keys())))
            for key in actual.keys():
                self.assertTrue(abs(actual[key] - expected[key]) < 0.02, "Unexpectedly large difference for %s in %s: %g versus %g" % (key, f_actual, actual[key], expected[key]))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
