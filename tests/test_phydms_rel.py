"""Runs the REL implementation of ``phydms``.

This test examines the functionality of the REL implementation of ``phydms``.

Written by Jonathan Mah
"""

import os
import unittest
import multiprocessing
import subprocess
import scipy
import pandas
import logging
import numpy
import glob


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

        outprefix = self.OUTPREFIX + self.MODEL + '_k2_' + str(self.K2)
        expected_prefix = self.EXPECTED_PREFIX + self.MODEL + '_k2_' + \
            str(self.K2)

        if self.MODEL is 'ExpCM':
            subprocess.check_call(['phydms', self.ALIGNMENT, self.TREE,
                                  'ExpCM_{0}'.format(self.PREFS), outprefix,
                                   '--ncpus', str(ncpus), '--gammaomega',
                                   '--empirical_bayes', str(self.K2)])
        elif self.MODEL is 'YNGKP_M5':
            subprocess.check_call(['phydms', self.ALIGNMENT, self.TREE,
                                  'YNGKP_M5', outprefix, '--ncpus',
                                   str(ncpus), '--empirical_bayes',
                                   str(self.K2)])
        else:
            raise ValueError('Only ExpCM and YNGKP models are implemented at '
                             'this time.')

        suffix_list = ['omegabycategory.csv', 'posteriorprobabilities.csv']
        for suffix in suffix_list:
            self.compare_output_dataframes(expected_prefix, outprefix, suffix)

        # Remove output files
        suffix_list.extend(['log.log', 'loglikelihood.txt', 'tree.newick',
                            'modelparams.txt'])
        for suffix in suffix_list:
            self.remove_output_files(outprefix, suffix)
        if os.path.isdir(self.OUTPREFIX):
            os.rmdir(self.OUTPREFIX)

    def compare_output_dataframes(self, expected_prefix, outprefix, suffix):
        values = {}
        file_list = [('expected', expected_prefix),
                     ('actual', outprefix)]
        for (name, prefix) in file_list:
            fname = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                    prefix + '_' + suffix))
            # fname = prefix + '_' + suffix
            self.assertTrue(
                os.path.isfile(fname),
                "Can't find output file {0}".format(fname))
            values[name] = pandas.read_csv(fname, comment='#', sep='\t')
        self.assertTrue(
            numpy.allclose(
                values['actual'].select_dtypes(exclude=[object]),
                values['expected'].select_dtypes(exclude=[object])) and
            values['actual'].select_dtypes(include=[object]).equals(
                values['expected'].select_dtypes(include=[object]))), \
            "Expected and actual results differ in value."

    def remove_output_files(self, outprefix, suffix):
        fname = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                outprefix + '_' + suffix))
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
