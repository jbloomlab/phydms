"""Runs ``phydms_comprehensive``

This test examines the entire functionality of ``phydms_comprehensive``
when run from the command line.

Written by Jesse Bloom and Sarah Hilton.
"""

import os
import unittest
import multiprocessing
import subprocess
import scipy
import pandas


class test_phydms_comprehensive(unittest.TestCase):
    """Tests command-line ``phydms_comprehensive``."""

    def test_NP(self):
        """Tests command-line ``phydms_comprehensive`` on NP data."""
        npdir = './NP_data/'
        prefs = '{0}/NP_prefs.tsv'.format(npdir)
        alignment = '{0}/NP_alignment.fasta'.format(npdir)
        tree = '{0}/NP_tree.newick'.format(npdir)
        for f in [prefs, alignment]:
            self.assertTrue(os.path.isfile(f), "Can't find file {0}".format(f))
        outprefix = './NP_test_results/'

        ncpus = min(20, multiprocessing.cpu_count())

        subprocess.check_call(['phydms_comprehensive', outprefix, alignment,
                prefs, "--tree", tree, "--omegabysite", '--brlen', 'scale',
                '--diffprefsbysite', '--ncpus', str(ncpus)])

        expectedresults = './expected_NP_test_results/'

        models = ['ExpCM_NP_prefs', 'averaged_ExpCM_NP_prefs', 'YNGKP_M0',
                'YNGKP_M5']

        for model in models:
            values = {}
            for (name, prefix) in [('expected', expectedresults), ('actual', outprefix)]:
                values[name] = {}
                for suffix in ['_loglikelihood.txt', '_modelparams.txt']:
                    with open(prefix + model + suffix) as f:
                        for line in f:
                            (x, y) = line.split('=')
                            values[name][x.strip()] = float(y)
            for param in values['actual'].keys():
                self.assertTrue(scipy.allclose(values['actual'][param],
                        values['expected'][param], atol=1e-2, rtol=1e-5))

            omegas = {}
            for (name, prefix) in [('expected', expectedresults), ('actual', outprefix)]:
                omegas[name] = pandas.read_csv(prefix + model + '_omegabysite.txt', 
                        comment='#', sep='\t')
                omegas[name] = omegas[name].sort_values(by='site', axis=0)
            self.assertTrue(scipy.allclose(omegas['actual']['P'].values,
                    omegas['expected']['P'].values, atol=0.01, rtol=0.03))
            sigsites = omegas['expected'][omegas['expected']['P'] < 0.05]['site'].values
            sigomegas = {}
            for (name, df) in omegas.items():
                sigomegas[name] = omegas[name][omegas[name]['site'].isin(sigsites)]['omega'].values
            self.assertTrue(((sigomegas['actual'] > 1) == (sigomegas['expected'] > 1)).all())

            if 'ExpCM' in model:
                diffprefs = {}
                largecutoff = 0.1
                smallcutoff = 0.01
                largesites = {}
                smallsites = {}
                for (name, prefix) in [('expected', expectedresults), ('actual', outprefix)]:
                    diffprefs[name] = pandas.read_csv(prefix + model + '_diffprefsbysite.txt',
                        comment='#', sep='\t')
                    # factor builds in a margin in expected versus actual 
                    factor = {'actual':1, 'expected':1.5}[name]
                    largesites[name] = set(diffprefs[name][diffprefs[name][
                            'half_sum_abs_dpi'] > largecutoff * factor]['site'].values)
                    smallsites[name] = set(diffprefs[name][diffprefs[name][
                            'half_sum_abs_dpi'] < smallcutoff / factor]['site'].values)
                self.assertTrue(largesites['actual'] >= largesites['expected'])
                self.assertTrue(smallsites['actual'] >= smallsites['expected'])


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
