"""Runs doctests.

Written by Jesse Bloom."""


import os
import sys
import glob
import doctest
import unittest
import doctest
import pkgutil
import phydmslib


def main():
    """Run the tests."""

    sys.stderr.write('Running tests...\n')
    failurestrings = []

    # test all modules with doctest
    for (importer, modname, ispkg) in pkgutil.iter_modules(phydmslib.__path__):
        if (not ispkg) and modname[0] != '_':
            sys.stderr.write('\nTesting %s with doctest... ' % modname)
            module = __import__('phydmslib.%s' % modname, None, None, modname.split('.'))
            suite = doctest.DocTestSuite(module)
            del module
            result = unittest.TestResult()
            suite.run(result)
            if result.wasSuccessful():
                sys.stderr.write('all %d tests were successful.\n' % result.testsRun)
            else:
                sys.stderr.write('test FAILED!\n')
                for (testcase, failstring) in result.failures:
                    failurestrings.append(failstring)

    # print summary of failures
    if not failurestrings:
        sys.stderr.write('\nTesting complete. All tests were passed successfully.\n')
    else:
        sys.stderr.write('\nTesting complete. Failed on the following tests:\n')
        for failstring in failurestrings:
            sys.stderr.write('\n*****************\n%s\n\n****************\n' % failstring)


if __name__ == '__main__':
    main() # run the program
