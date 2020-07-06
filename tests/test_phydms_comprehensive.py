"""Runs ``phydms_comprehensive``

This test examines the functionality of ``phydms_comprehensive`` when run
from the command-line with the `--gammaomega` and
`--omega_random_effects_likelihood` flags.

Written by Jesse Bloom, Sarah Hilton, and Jonathan Mah
"""

import os
import unittest
import multiprocessing
import subprocess
import scipy
import pandas


class test_phydms_comprehensive(unittest.TestCase):
    """Tests command-line ``phydms_comprehensive`` with the `--gammaomega`
    and `--omega_random_effects_likelihood` flags. This test is performed with
    the minimum number of categories used for integration, being 2."""

    def test_NP(self):
        """Tests command-line ``phydms_comprehensive`` on NP data."""
        tree = os.path.abspath(os.path.join(
            os.path.dirname(__file__), './NP_data/NP_tree_short.newick'))
        alignment = os.path.abspath(os.path.join(
            os.path.dirname(__file__), './NP_data/NP_alignment_short.fasta'))
        prefs = os.path.abspath(os.path.join(
            os.path.dirname(__file__), './NP_data/NP_prefs_short.csv'))
        for f in [prefs, alignment]:
            self.assertTrue(os.path.isfile(f), "Can't find file {0}".format(f))
        outprefix = os.path.abspath(os.path.join(
            os.path.dirname(__file__), './NP_test_results/'))
        if outprefix[-1] != "/":
            outprefix = "{0}/".format(outprefix)

        ncpus = min(20, multiprocessing.cpu_count())

        K = 2  # Number of bins used in fitting gamma distribution

        J = 2  # Number of bins used in empirical Bayesian integration

        subprocess.check_call(
            ['phydms_comprehensive', outprefix, alignment,
             prefs, "--tree", tree, "--omegabysite", '--brlen', 'scale',
             '--ncpus', str(ncpus), '--gammaomega', '--ncats', str(K),
             '--omega_random_effects_likelihood', '--REL_ncats', str(J)])

        expectedresults = os.path.abspath(os.path.join(
            os.path.dirname(__file__), './expected_NP_test_results/'))

        models = ['ExpCM_NP_prefs_short', 'YNGKP_M0', 'YNGKP_M5',
                  'ExpCM_NP_prefs_short_gammaomega']

        for model in models:
            values = {}
            for (name, prefix) in [('expected', expectedresults),
                                   ('actual', outprefix)]:
                values[name] = {}
                for suffix in ['_loglikelihood.txt', '_modelparams.txt']:
                    fname = os.path.abspath(os.path.join(
                        prefix, './{0}{1}'.format(model, suffix)))
                    with open(fname) as f:
                        for line in f:
                            (x, y) = line.split('=')
                            values[name][x.strip()] = float(y)
            for param in values['actual'].keys():
                self.assertTrue(scipy.allclose(
                    values['actual'][param],
                    values['expected'][param], atol=1e-2, rtol=1e-5))
            # FEL
            omegas = {}
            for (name, prefix) in [('expected', expectedresults),
                                   ('actual', outprefix)]:
                fname = os.path.abspath(os.path.join(
                    prefix, './{0}{1}'.format(model, '_omegabysite.txt')))
                omegas[name] = pandas.read_csv(
                    fname, comment='#', sep='\t')
                omegas[name] = omegas[name].sort_values(by='site', axis=0)
            self.assertTrue(scipy.allclose(
                omegas['actual']['P'].values,
                omegas['expected']['P'].values, atol=0.01, rtol=0.03))
            sigsites = omegas['expected'][
                omegas['expected']['P'] < 0.05]['site'].values
            sigomegas = {}
            for (name, df) in omegas.items():
                sigomegas[name] = omegas[name][
                    omegas[name]['site'].isin(sigsites)]['omega'].values
            self.assertTrue(((
                sigomegas['actual'] > 1) == (sigomegas['expected'] > 1)).all())

        # REL
            if 'gammomega' in model or model == 'YNGKP_M5':
                omegas = {}
                for (name, prefix) in [('expected', expectedresults),
                                       ('actual', outprefix)]:
                    fname = os.path.abspath(os.path.join(
                        prefix, './{0}{1}'.format(model, '_omegabycategory.csv')))
                    omegas[name] = pandas.read_csv(fname)
                self.assertTrue(scipy.allclose(
                    omegas['actual']['post_probability'].values,
                    omegas['expected']['post_probability'].values,
                    atol=0.001, rtol=0.003))
                self.assertTrue(scipy.allclose(
                    omegas['actual']['omega'].values,
                    omegas['expected']['omega'].values,
                    atol=0.001, rtol=0.003))

                posteriors = {}
                for (name, prefix) in [('expected', expectedresults),
                                       ('actual', outprefix)]:
                    fname = os.path.abspath(os.path.join(
                        prefix, './{0}{1}'.format(model, '_posteriorprobabilities.csv')))
                    posteriors[name] = pandas.read_csv(fname)
                self.assertTrue(scipy.allclose(
                    posteriors['actual']['p(omega > 1)'].values,
                    posteriors['expected']['p(omega > 1)'].values,
                    atol=0.001, rtol=0.003))
                self.assertTrue(scipy.allclose(
                    posteriors['actual']['fdr'].values,
                    posteriors['expected']['fdr'].values,
                    atol=0.001, rtol=0.003))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
