"""Tests ``phydms_prepalignment``.

Written by Jesse Bloom."""


import sys
import os
import unittest
import subprocess


class TestPrepAlignment(unittest.TestCase):
    """
    Runs ``phydms_prepalignment`` on test data, compares to known results.
    """

    def setUp(self):
        """Gets files set up appropriately."""
        self.testdir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                './prepalignment_tests/'))


    def test_unaligned(self):
        """Test on data set that needs aligning."""
        inseqs = os.path.abspath(os.path.join(self.testdir,
                'unaligned_SUMO1_orthologs.fasta'))
        for (minidentity, minuniqueness, purgeseqs, keepseqs) in [
                ('0.7', '2', [''], [''])
                ]:
            alignmentSuffix = 'test_unaligned_SUMO1_alignment_minidentity{0}_minuniqueness{1}_purgeseqs{2}_keepseqs{3}.fasta'.format(minidentity, minuniqueness, os.path.splitext(os.path.basename(''.join(purgeseqs)))[0], os.path.splitext(os.path.basename(''.join(keepseqs)))[0])
            alignment = os.path.abspath(os.path.join(self.testdir, alignmentSuffix))
            subprocess.check_call([
                    'phydms_prepalignment',
                    inseqs,
                    alignment,
                    'Hsap',
                    '--minidentity', minidentity,
                    '--minuniqueness', minuniqueness,
                    '--purgeseqs'] + purgeseqs +
                    ['--keepseqs'] + keepseqs
                )
            expected = alignment.replace('test_SUMO', 'expected_SUMO')
            with open(expected) as f_expected, open(alignment) as f_actual:
                expectedtext = f_expected.read()
                actualtext = f_actual.read()
            self.assertTrue(expectedtext == actualtext, "Did not get expected results for {0}".format(alignment))


    def test_prealigned(self):
        """Test on pre-aligned data set."""
        inseqs = os.path.abspath(os.path.join(self.testdir,
                'prealigned_SUMO1_orthologs.fasta'))
        purgeseqsfile = os.path.abspath(os.path.join(self.testdir,
                'test_purgeseqsfile.txt'))
        with open(purgeseqsfile, 'w') as f:
            f.write('Fcat\n\n')
        keepseqsfile = os.path.abspath(os.path.join(self.testdir,
                'test_keepseqsfile.txt'))
        with open(keepseqsfile, 'w') as f:
            f.write('Panu\nGgor\n')
        for (minidentity, minuniqueness, purgeseqs, keepseqs) in [
                ('0.7', '2', [''], ['']),
                ('0.6', '2', [''], ['']),
                ('0.7', '1', [''], ['']),
                ('0.7', '2', [''], ['Panu', 'Ggor']),
                ('0.7', '2', [''], [keepseqsfile]),
                ('0.7', '2', [purgeseqsfile], ['Panu', 'Ggor']),
                ('0.7', '2', ['Fcat'], [keepseqsfile]),
                ]:
            alignmentSuffix = 'test_SUMO1_alignment_minidentity{0}_minuniqueness{1}_purgeseqs{2}_keepseqs{3}.fasta'.format(minidentity, minuniqueness, os.path.splitext(os.path.basename(''.join(purgeseqs)))[0], os.path.splitext(os.path.basename(''.join(keepseqs)))[0])
            alignment = os.path.abspath(os.path.join(self.testdir, alignmentSuffix))
            subprocess.check_call([
                    'phydms_prepalignment',
                    inseqs,
                    alignment,
                    'Hsap',
                    '--prealigned',
                    '--minidentity', minidentity,
                    '--minuniqueness', minuniqueness,
                    '--purgeseqs'] + purgeseqs +
                    ['--keepseqs'] + keepseqs
                )
            expected = alignment.replace('test_SUMO', 'expected_SUMO')
            with open(expected) as f_expected, open(alignment) as f_actual:
                expectedtext = f_expected.read()
                actualtext = f_actual.read()
            self.assertTrue(expectedtext == actualtext, "Did not get expected results for {0}".format(alignment))


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    unittest.main(testRunner=runner)
