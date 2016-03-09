"""Tests ``phydms`` on a dataset of Gal4 sequences.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess
import dms_tools.file_io



class TestOnGal4(unittest.TestCase):
    """
    Runs ``phydms_comprehensive`` on test data, compares to known results.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.expected_dir = './expected_Gal4_results/'
        self.test_dir = './Gal4_test_results/'
        self.prefs = 'Gal4_prefs.txt'
        self.assertTrue(os.path.isfile(self.prefs), "Can't find required file %s" % self.prefs)
        self.alignment = 'Gal4s.fasta'
        self.assertTrue(os.path.isfile(self.alignment), "Can't find required file %s" % self.alignment)
        self.all_models = ['ExpCM_Gal4_prefs', 'averaged_ExpCM_Gal4_prefs', 'YNGKP_M0', 'YNGKP_M3']
        self.likelihood_files = ['%s_loglikelihood.txt' % model for model in self.all_models]
        self.params_files = ['%s_modelparams.txt' % model for model in self.all_models]
        self.bysite_files =\
                ['%s_omegabysite.txt' % model for model in self.all_models if 'YNGKP_M0' not in model]
        self.diffpref_files = ['%s_diffprefsbysite.txt' % model for model in self.all_models if 'YNGKP' not in model]
        for f in self.likelihood_files + self.params_files + self.bysite_files + self.diffpref_files:
            fname = "%s/%s" % (self.expected_dir, f)
            self.assertTrue(os.path.isfile(fname), "Cannot find required file %s" % fname)
            toremove = "%s/%s" % (self.test_dir, f)
            if os.path.isfile(toremove):
                os.remove(toremove)


    def test_Phydms(self):
        """Runs ``phydms_comprehensive``."""
        cmds = ['phydms_comprehensive', self.test_dir, self.alignment, self.prefs, '--ncpus', '-1', '--yngkp', 'M3', '--no-stringencybysite', '--fitF3X4']
        sys.stderr.write('\nRunning phydms with the the following command:\n%s\n' % ' '.join(cmds))
        subprocess.call(cmds)

        sys.stderr.write('\nTesting for presence of expected output files...\n')
        for f in self.likelihood_files + self.params_files + self.bysite_files + self.diffpref_files:
            fname = '%s/%s' % (self.test_dir, f)
            self.assertTrue(os.path.isfile(fname), "Failed to created expected file %s" % fname)

        sys.stderr.write('\nTesting for expected likelihood values...\n')
        for f in self.likelihood_files:
            with open('%s/%s' % (self.expected_dir, f)) as fin:
                expected = float(fin.read().split('=')[-1])
            with open('%s/%s' % (self.test_dir, f)) as fin:
                actual = float(fin.read().split('=')[-1])
            self.assertTrue(abs(expected - actual) < 1, "Unexpectedly large differences in %s: %g versus %g" % (f, expected, actual))

        sys.stderr.write('\nTesting for expected parameter values...\n')
        for f in self.params_files:
            with open('%s/%s' % (self.expected_dir, f)) as fin:
                expected = dict([(line.split('=')[0], float(line.split('=')[1])) for line in fin.readlines()])
            with open('%s/%s' % (self.test_dir, f)) as fin:
                actual = dict([(line.split('=')[0], float(line.split('=')[1])) for line in fin.readlines()])
            self.assertTrue(set(expected.keys()) == set(actual.keys()), "Different parameters in %s: %s versus %s" % (f, str(expected.keys()), str(actual.keys())))
            for key in actual.keys():
                self.assertTrue(abs(actual[key] - expected[key]) < 0.02, "Unexpectedly large difference for %s in %s: %g versus %g" % (key, f, actual[key], expected[key]))

        sys.stderr.write('\nTesting for expected selection by site values...\n')
        for f in self.bysite_files:
            with open('%s/%s' % (self.expected_dir, f)) as fin:
                expected = dict([(line.split()[0], {'P':float(line.split()[2]), 'value':float(line.split()[1])}) for line in fin.readlines() if line[0] != '#'])
            with open('%s/%s' % (self.test_dir, f)) as fin:
                actual = dict([(line.split()[0], {'P':float(line.split()[2]), 'value':float(line.split()[1])}) for line in fin.readlines() if line[0] != '#'])
            self.assertTrue(set(expected.keys()) == set(actual.keys()), "Different sites in %s" % f)
            for site in actual.keys():
                self.assertTrue(abs(actual[site]['P'] - expected[site]['P']) < 0.025 * max(actual[site]['P'], expected[site]['P']), "Unexpectedly large difference in P for %s in %s: %g versus %g" % (site, f, actual[site]['P'], expected[site]['P']))
                if actual[site]['P'] < 0.1:
                    self.assertTrue(abs(actual[site]['value'] - expected[site]['value']) < 0.025 * max(actual[site]['value'], expected[site]['value']), "Unexpectedly large difference in value for %s in %s: %g versus %g" % (site, f, actual[site]['value'], expected[site]['value']))

        sys.stderr.write("\nTesting for expected differential preferences at each site...\n")
        for f in self.diffpref_files:
            actual = dms_tools.file_io.ReadDiffPrefs('%s/%s' % (self.test_dir, f))[2]
            expected = dms_tools.file_io.ReadDiffPrefs('%s/%s' % (self.expected_dir, f))[2]
            self.assertTrue(set(expected.keys()) == set(actual.keys()), "Different sites in %s" % f)
            for site in actual.keys():
                for aa in actual[site].keys():
                    (x, y) = (actual[site][aa], expected[site][aa])
                    self.assertTrue(abs(x - y) < 0.02, "Unexpectedly large difference in diffpref for site %s, amino-acid %s in %s: got %g, expected %g" % (site, aa, f, x, y))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
