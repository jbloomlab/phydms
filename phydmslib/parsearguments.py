"""
Module for parsing arguments.
"""


import sys
import os
import re
import argparse
import phydmslib


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""
    def error(self, message):
        """Prints error message, then help."""
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def NonNegativeInt(n):
    """If *n* is non-negative integer returns it, otherwise an error.

    >>> print "%d" % NonNegativeInt('8')
    8

    >>> NonNegativeInt('8.1')
    Traceback (most recent call last):
       ...
    ArgumentTypeError: 8.1 is not an integer

    >>> print "%d" % NonNegativeInt('0')
    0

    >>> NonNegativeInt('-1')
    Traceback (most recent call last):
       ...
    ArgumentTypeError: -1 is not non-negative

    """
    if not isinstance(n, str):
        raise argparse.ArgumentTypeError('%r is not a string' % n)
    try:
       n = int(n)
    except:
        raise argparse.ArgumentTypeError('%s is not an integer' % n)
    if n < 0:
        raise argparse.ArgumentTypeError('%d is not non-negative' % n)
    else:
        return n


def FloatGreaterThanZero(x):
    """If *x* is string for float > 0, returns it, otherwise an error.

    Designed based on this: http://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin

    >>> print "%.3f" % FloatGreaterThanZero('0.1')
    0.100

    >>> FloatGreaterThanZero('0.0')
    Traceback (most recent call last):
        ...
    ArgumentTypeError: 0.0 not a float greater than zero

    >>> FloatGreaterThanZero('hi')
    Traceback (most recent call last):
        ...
    ValueError: could not convert string to float: hi
    """
    x = float(x)
    if x > 0:
        return x
    else:
        raise argparse.ArgumentTypeError("%r not a float greater than zero" % x)


def ExistingFile(fname):
    """If *fname* is name of an existing file return it, otherwise an error."""
    if os.path.isfile(fname):
        return fname
    else:
        raise argparse.ArgumentTypeError("%s does not specify a valid file name" % fname)


def TreeFile(fname):
    """Returns *fname* if an existing file or *random* or *nj*, otherwise an error."""
    if os.path.isfile(fname):
        if fname.lower() in ['random', 'nj']:
            raise argparse.ArgumentTypeError("Ambiguous meaning of tree since there is an existing file named %s" % fname)
        else:
            return fname
    elif fname.lower() in ['random', 'nj']:
        return fname.lower()
    else:
        raise argparse.ArgumentTypeError("Invalid value for tree: must be existing file, 'nj', or 'random'.")


def ModelOption(model):
    """Returns *model* if a valid choice.
    
    Returns the string if it specifies a ``YNGKP_`` model variant.
    
    Returns *('ExpCM', prefsfile)* if it specifies an ``ExpCM_`` model.
    """
    if model in ['YNGKP_M0', 'YNGKP_M7', 'YNGKP_M8']:
        return model
    elif len(model) > 6 and model[ : 6] == 'ExpCM_':
        fname = model[6 : ] 
        if os.path.isfile(fname):
            return ('ExpCM', fname)
        else:
            raise argparse.ArgumentTypeError("ExpCM_ must be followed by the name of an existing file. You specified the following, which is not an existing file: %s" % fname)
    else:
        raise argparse.ArgumentTypeError("Invalid model")


def PhyDMSParser():
    """Returns *argparse.ArgumentParser* for ``phydms`` script."""
    parser = ArgumentParserNoArgHelp(description='Phylogenetic inference using deep mutational scanning data. Version %s by %s. Full documentation at %s' % (phydmslib.__version__, phydmslib.__author__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('alignment', help='Existing FASTA file with aligned codon sequences.', type=ExistingFile)
    parser.add_argument('tree', help="Existing Newick tree file or 'random' or 'nj'.", type=TreeFile)
    parser.add_argument('model', help="Codon substitution model: YNGKP_M0, YNGKP_M7, YNGKP_M8, ExpCM_<prefsfilename>", type=ModelOption)
    parser.add_argument('outprefix', help='Prefix for output files.', type=str)
    parser.set_defaults(omegabysite=False)
    parser.add_argument('--omegabysite', dest='omegabysite', action='store_true', help="Fit a different omega (dN/dS) for each site?")
    parser.set_defaults(infertopology=False)
    parser.add_argument('--infertopology', dest='infertopology', action='store_true', help="Infer topology starting from 'tree'; otherwise topology fixed to 'tree'. Requires '--oldlikmethod'.")
    parser.set_defaults(oldlikmethod=False)
    parser.add_argument('--oldlikmethod', dest='oldlikmethod', action='store_true', help='Use old Bio++ likelihood method. Only allowed for non-partitioned models.')
    parser.set_defaults(fixbrlen=False)
    parser.add_argument('--fixbrlen', dest='fixbrlen', action='store_true', help="Fix branch lengths to those of initial 'tree'?")
    parser.set_defaults(fitF3X4=False)
    parser.add_argument('--fitF3X4', dest='fitF3X4', action='store_true', help='Fit F3X4 frequencies for YNGKP; otherwise use empirical')
    parser.add_argument('--minbrlen', default=1e-6, help='Min branch length for starting tree.', type=FloatGreaterThanZero)
    parser.add_argument('--seed', type=int, default=1, help="Random number seed.")
    parser.set_defaults(recursion='S')
    # comment out this option for now as the double recursion seems not to work properly
    #parser.add_argument('--recursion', choices=['S', 'D'], default='S', help='Likelihood recursion for Bio++.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=phydmslib.__version__))
    return parser


if __name__ == '__main__':
    import doctest
    doctest.testmod()
