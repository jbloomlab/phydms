"""Runs the REL implementation of ``phydms``.

This test examines the functionality of the REL implementation of ``phydms``.

Written by Jonathan Mah
"""

import os
import unittest
import multiprocessing
import subprocess
import pandas
import scipy


class test_phydms_rel_ExpCM_k2_4(unittest.TestCase):
    """Tests command-line usage of REL implementation of ``phydms`` with 4 bins
       used for integration.
    """
    TREE = os.path.abspath(os.path.join(os.path.dirname(__file__),
                           './REL_input_data/NP_tree.newick'))
    ALIGNMENT = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                './REL_input_data/NP_alignment.fasta'))
    PREFS = os.path.abspath(os.path.join(os.path.dirname(__file__),
                            './REL_input_data/NP_prefs.csv'))
    MODEL = 'ExpCM'
    K2 = 4  # Number of bins used for empirical_bayes integration
    OUTPREFIX = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                './_phydms_rel_test_results/'))
    EXPECTED_PREFIX = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                      './expected_rel_test_results/'))

    if OUTPREFIX[-1] != '/':
        OUTPREFIX = '{0}/'.format(OUTPREFIX)
    if EXPECTED_PREFIX[-1] != '/':
        EXPECTED_PREFIX = '{0}/'.format(EXPECTED_PREFIX)

    def test_rel_expected_output_ExpCM(self):
        """Tests expected output from command-line usage of REL implementation
        of ``phydms`` using NP data for ExpCM.
        """
        # Confirm input data exist
        for f in [self.ALIGNMENT, self.TREE, self.PREFS]:
            self.assertTrue(
                os.path.isfile(f), "Can't find test file {0}".format(f))

        ncpus = 1
        model_with_bins = self.MODEL + '_k2_' + str(self.K2)

        if self.MODEL is 'ExpCM':
            subprocess.check_call(['phydms', self.ALIGNMENT, self.TREE,
                                  'ExpCM_{0}'.format(self.PREFS),
                                   self.OUTPREFIX + model_with_bins,
                                   '--ncpus', str(ncpus), '--gammaomega',
                                   '--empirical_bayes', str(self.K2)])
        elif self.MODEL is 'YNGKP_M5':
            subprocess.check_call(['phydms', self.ALIGNMENT, self.TREE,
                                  'YNGKP_M5', self.OUTPREFIX + model_with_bins,
                                   '--ncpus', str(ncpus), '--empirical_bayes',
                                   str(self.K2)])
        else:
            raise ValueError('Only ExpCM and YNGKP models are implemented at '
                             'this time.')

        self.compare_output_dataframes(
            self.OUTPREFIX, self.EXPECTED_PREFIX, model_with_bins)

        self.remove_output_files(self.OUTPREFIX, model_with_bins)
        if not os.listdir(self.OUTPREFIX):
            os.rmdir(self.OUTPREFIX)

    def compare_output_dataframes(self, outprefix, expected_prefix,
                                  model_with_bins):
        omegas = {}
        for (name, prefix) in [('expected', expected_prefix),
                               ('actual', outprefix)]:
            fname = os.path.abspath(os.path.join(prefix, './{0}{1}'.format(
                model_with_bins, '_omegabycategory.csv')))
            omegas[name] = pandas.read_csv(fname)
        self.assertTrue(scipy.allclose(
            omegas['actual']['post_probability'].values,
            omegas['expected']['post_probability'].values,
            atol=0.001, rtol=0.003))
        self.assertTrue(scipy.allclose(
            omegas['actual']['omega_value'].values,
            omegas['expected']['omega_value'].values,
            atol=0.001, rtol=0.003))

        posteriors = {}
        for (name, prefix) in [('expected', expected_prefix),
                               ('actual', outprefix)]:
            fname = os.path.abspath(os.path.join(prefix, './{0}{1}'.format(
                model_with_bins, '_posteriorprobabilities.csv')))
            posteriors[name] = pandas.read_csv(fname)
        self.assertTrue(scipy.allclose(
            posteriors['actual']['pr(positive_selection)'].values,
            posteriors['expected']['pr(positive_selection)'].values,
            atol=0.001, rtol=0.003))

    def remove_output_files(self, outprefix, model_with_bins):
        suffix_list = ['_omegabycategory.csv', '_posteriorprobabilities.csv',
                       '_log.log', '_loglikelihood.txt', '_tree.newick',
                       '_modelparams.txt']
        for suffix in suffix_list:
            fname = os.path.abspath(os.path.join(outprefix, './{0}{1}'.format(
                model_with_bins, suffix)))
            if os.path.isfile(fname):
                os.remove(fname)


class test_phydms_rel_ExpCM_k2_50(test_phydms_rel_ExpCM_k2_4):
    """Tests command-line usage of REL implementation of ``phydms`` with 50
       bins used for integration.
    """
    K2 = 50


class test_phydms_rel_YNGKP_M5_k2_4(test_phydms_rel_ExpCM_k2_4):
    """Tests command-line usage of REL implementation of YNGKP_M5 with 4
       bins used for integration.
    """
    MODEL = 'YNGKP_M5'


class test_phydms_rel_YNGKP_M5_k2_50(test_phydms_rel_YNGKP_M5_k2_4):
    """Tests command-line usage of REL implementation of YNGKP_M5 with 50
       bins used for integration.
    """
    K2 = 50


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
