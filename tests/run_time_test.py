"""Runs all tests."""


import os
import sys
import glob
import unittest
import pandas as pd
import time

def main():
    """Runs the tests."""
    r = {"test":[], "time":[]}
    failurestrings = []
    for test in glob.glob('test_*.py'):
        s = time.time()
        sys.stderr.write('\nRunning tests in {0}...\n'.format(test))
        test = os.path.splitext(test)[0]
        suite = unittest.TestLoader().loadTestsFromName(test)
        result = unittest.TestResult()
        suite.run(result)
        e = time.time()
        if result.wasSuccessful():
            sys.stderr.write('All tests were successful.\n')
            r["test"].append(test)
            r["time"].append(float(e-s)/60)
        else:
            sys.stderr.write('Test(s) FAILED!\n')
            for (testcase, failstring) in result.failures + result.errors:
                failurestrings.append(failstring)

    if not failurestrings:
        sys.stderr.write('\nTesting complete. All passed successfully.\n')
    else:
        sys.stderr.write('\nTesting complete. Failed on the following:\n')
        for fstring in failurestrings:
            sys.stderr.write('\n*********\n{0}\n********\n'.format(fstring))
    r = pd.DataFrame(r)
    r["Greater_10min"] = ["y" if x>10 else "n" for x in r["time"]]
    r.to_csv("time_test_results.csv", index=False)

if __name__ == '__main__':
    main()
