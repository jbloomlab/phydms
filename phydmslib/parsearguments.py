"""
Module for parsing arguments.
"""


import sys
import os
import re
import argparse
import phydmslib
import phydmslib.constants


# allowed variants of YNGKP models
yngkp_modelvariants = ['M0', 'M5']

class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""
    def error(self, message):
        """Prints error message, then help."""
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)

class ArgumentDefaultsRawDescriptionFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    """Print default arguments and raw text description formatter.

    Based on this:
    http://stackoverflow.com/questions/18462610/argumentparser-epilog-and-description-formatting-in-conjunction-with-argumentdef
    """
    pass


def NonNegativeInt(n):
    """If *n* is non-negative integer returns it, otherwise an error.

    >>> print("%d" % NonNegativeInt('8'))
    8

    >>> NonNegativeInt('8.1')
    Traceback (most recent call last):
       ...
    ValueError: 8.1 is not an integer

    >>> print("%d" % NonNegativeInt('0'))
    0

    >>> NonNegativeInt('-1')
    Traceback (most recent call last):
       ...
    ValueError: -1 is not non-negative

    """
    if not isinstance(n, str):
        raise ValueError('%r is not a string' % n)
    try:
       n = int(n)
    except:
        raise ValueError('%s is not an integer' % n)
    if n < 0:
        raise ValueError('%d is not non-negative' % n)
    else:
        return n

def IntGreaterThanZero(n):
    """If *n* is an integer > 0, returns it, otherwise an error."""
    try:
        n = int(n)
    except:
        raise ValueError("%s is not an integer" % n)
    if n <= 0:
        raise ValueError("%d is not > 0" % n)
    else:
        return n

def IntGreaterThanOne(n):
    """If *n* is an integer > 1, returns it, otherwise an error."""
    try:
        n = int(n)
    except:
        raise ValueError("%s is not an integer" % n)
    if n <= 1:
        raise ValueError("%d is not > 1" % n)
    else:
        return n

def FloatGreaterThanEqualToZero(x):
    """If *x* is a float >= 0, returns it, otherwise raises and error.

    >>> print('%.1f' % FloatGreaterThanEqualToZero('1.5'))
    1.5

    >>> print('%.1f' % FloatGreaterThanEqualToZero('-1.1'))
    Traceback (most recent call last):
       ...
    ValueError: -1.1 not float greater than or equal to zero
    """
    try:
        x = float(x)
    except:
        raise ValueError("%r not float greater than or equal to zero" % x)
    if x >= 0:
        return x
    else:
        raise ValueError("%r not float greater than or equal to zero" % x)


def FloatGreaterThanOne(x):
    """If *x* is a string for a float > 1, returns it, otherwise an error."""
    x = float(x)
    if x > 1:
        return x
    else:
        raise ValueError("%r not a float greater than one" % x)


def FloatGreaterThanZero(x):
    """If *x* is string for float > 0, returns it, otherwise an error.

    Designed based on this: http://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin

    >>> print("%.3f" % FloatGreaterThanZero('0.1'))
    0.100

    >>> FloatGreaterThanZero('0.0')
    Traceback (most recent call last):
        ...
    ValueError: 0.0 not a float greater than zero
    """
    x = float(x)
    if x > 0:
        return x
    else:
        raise ValueError("%r not a float greater than zero" % x)


def FloatBetweenZeroAndOne(x):
    """Returns *x* only if *0 <= x <= 1*, otherwise raises error."""
    x = float(x)
    if 0 <= x <= 1:
        return x
    else:
        raise ValueError("{0} not a float between 0 and 1.".format(x))


def diffPrefsPrior(priorstring):
    """Parses `priorstring` and returns `prior` tuple."""
    assert isinstance(priorstring, str)
    prior = priorstring.split(',')
    if len(prior) == 3 and prior[0] == 'invquadratic':
        [c1, c2] = [float(x) for x in prior[1 : ]]
        assert c1 > 0 and c2 > 0, "C1 and C2 must be > 1 for invquadratic prior"
        return ('invquadratic', c1, c2)
    else:
        raise ValueError("Invalid diffprefsprior: {0}".format(priorstring))


def ExistingFile(fname):
    """If *fname* is name of an existing file return it, otherwise an error."""
    if os.path.isfile(fname):
        return fname
    else:
        raise ValueError("%s must specify a valid file name" % fname)

def ExistingFileOrNone(fname):
    """Like `Existingfile`, but if `fname` is string "None" then return `None`."""
    if os.path.isfile(fname):
        return fname
    elif fname.lower() == 'none':
        return None
    else:
        raise ValueError("%s must specify a valid file name or 'None'" % fname)


def ModelOption(model):
    """Returns *model* if a valid choice.

    Returns the string if it specifies a ``YNGKP_`` model variant.

    Returns *('ExpCM', prefsfile)* if it specifies an ``ExpCM_`` model.
    """
    yngkpmatch = re.compile('^YNGKP_M[{0}]$'.format(''.join([m[1 : ] for m in yngkp_modelvariants])))
    if yngkpmatch.search(model):
        return model
    elif len(model) > 6 and model[ : 6] == 'ExpCM_':
        fname = model[6 : ]
        if os.path.isfile(fname):
            return ('ExpCM', fname)
        else:
            raise ValueError("ExpCM_ must be followed by the name of an existing file. You specified the following, which is not an existing file: %s" % fname)
    else:
        raise ValueError("Invalid model")


def PhyDMSPrepAlignmentParser():
    """Returns *argparse.ArgumentParser* for ``phydms_prepalignment``."""
    parser = ArgumentParserNoArgHelp(formatter_class=ArgumentDefaultsRawDescriptionFormatter,
            description='\n'.join([
            "Prepare alignment of protein-coding DNA sequences.\n",
            "Steps:",
            " * Any sequences specified by '--purgeseqs' are removed.",
            " * Sequences not of length divisible by 3 are removed.",
            " * Sequences with ambiguous nucleotides are removed.",
            " * Sequences with non-terminal stop codons are removed;",
            "   terminal stop codons are trimmed.",
            " * Sequences that do not encode unique proteins are removed",
            "   unless they are specified for retention by '--keepseqs'.",
            " * A multiple sequence alignment is built using MAFFT.",
            "   This step is skipped if you specify '--prealigned'.",
            " * Sites gapped in reference sequence are stripped.",
            " * Sequences with too little protein identity to reference",
            "   sequence are removed, counting both mismatches and unstripped",
            "   gaps as differences. Identity cutoff set by '--minidentity'.",
            " * Sequences too similar to other sequences are removed. An",
            "   effort is made to keep one representative of sequences found",
            "   many times in input set. Uniqueness threshold set ",
            "   by '--minuniqueness'. You can specify sequences to not",
            "   remove via '--keepseqs'.",
            " * Problematic characters in header names are replaced by",
            "   underscores. This is any space, comma, colon, semicolon",
            "   parenthesis, bracket, single quote, or double quote.",
            " * An alignment is written, as well as a plot with same root",
            "   but extension '.pdf' that shows divergence from reference",
            "   of all sequences retained and purged due to identity or",
            "   uniqueness.\n",
            phydmslib.__acknowledgments__,
            'Version {0}'.format(phydmslib.__version__),
            'Full documentation at {0}'.format(phydmslib.__url__),
            ]))
    parser.add_argument('inseqs', type=ExistingFile, help="FASTA file giving input coding sequences.")
    parser.add_argument('alignment', help="Name of created output FASTA alignment. PDF plot has same root, but extension '.pdf'.")
    parser.add_argument('refseq', help="Reference sequence in 'inseqs': specify substring found ONLY in header for that sequence.")
    parser.set_defaults(prealigned=False)
    parser.add_argument('--prealigned', action='store_true', dest='prealigned', help="Sequences in 'inseqs' are already aligned, do NOT re-align.")
    parser.add_argument('--mafft', help="Path to MAFFT (http://mafft.cbrc.jp/alignment/software/).", default='mafft')
    parser.add_argument('--minidentity', type=FloatBetweenZeroAndOne, help="Purge sequences with <= this protein identity to 'refseq'.", default=0.7)
    parser.add_argument('--minuniqueness', type=IntGreaterThanZero, default=2, help="Require each sequence to have >= this many protein differences relative to other sequences.")
    parser.add_argument('--purgeseqs', nargs='*', help="Specify sequences to always purge. Any sequences with any of the substrings specified here are always removed. The substrings can either be passed as repeated arguments here, or as the name of an existing file which has one substring per line.")
    parser.add_argument('--keepseqs', nargs='*', help="Do not purge any of these sequences for lack of identity or uniqueness. Specified in the same fashion as for '--purgeseqs'.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=phydmslib.__version__))
    return parser


def PhyDMSLogoPlotParser():
    """Returns `argparse.ArgumentParser` for ``phydms_logoplot``."""
    parser = ArgumentParserNoArgHelp(description=
            "Make logo plot of preferences or differential preferences. "
            "Uses weblogo (http://weblogo.threeplusone.com/). "
            "{0} Version {1}. Full documentation at {2}".format(
            phydmslib.__acknowledgments__,
            phydmslib.__version__, phydmslib.__url__),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prefs', type=ExistingFile, help="File with "
            "amino-acid preferences; same format as input to 'phydms'.")
    group.add_argument('--diffprefs', type=ExistingFile, help="File with "
            "differential preferences; in format output by 'phydms'.")
    parser.add_argument('outfile', help='Name of created PDF logo plot.')
    parser.add_argument('--stringency', type=FloatGreaterThanEqualToZero,
            default=1, help="Stringency parameter to re-scale prefs.")
    parser.add_argument('--nperline', type=IntGreaterThanZero, default=70,
            help="Number of sites per line.")
    parser.add_argument('--numberevery', type=IntGreaterThanZero, default=10,
            help="Number sites at this interval.")
    parser.add_argument('--mapmetric', default='functionalgroup', choices=['kd',
            'mw', 'charge', 'functionalgroup'], help='Metric used to color '
            'amino-acid letters. kd = Kyte-Doolittle hydrophobicity; '
            'mw = molecular weight; functionalgroup = divide in 7 '
            'groups; charge = charge at neutral pH.')
    parser.add_argument('--colormap', type=str, default='jet',
            help="Name of `matplotlib` color map for amino acids "
            "when `--mapmetric` is 'kd' or 'mw'.")
    parser.add_argument('--diffprefheight', type=FloatGreaterThanZero,
            default=1.0, help="Height of diffpref logo in each direction.")
    parser.add_argument('--omegabysite', help="Overlay omega on "
            "logo plot. Specify '*_omegabysite.txt' file from 'phydms'.",
            type=ExistingFileOrNone)
    parser.add_argument('--minP', type=FloatGreaterThanZero, default=1e-4,
            help="Min plotted P-value for '--omegabysite' overlay.")
    parser.add_argument('-v', '--version', action='version',
            version='%(prog)s {version}'.format(version=phydmslib.__version__))
    return parser


def PhyDMSComprehensiveParser():
    """Returns *argparse.ArgumentParser* for ``phdyms_comprehensive`` script."""
    parser = ArgumentParserNoArgHelp(description=("Comprehensive phylogenetic "
            "model comparison and detection of selection informed by deep "
            "mutational scanning data. This program runs 'phydms' repeatedly "
            "to compare substitution models and detect selection. "
            "{0} Version {1}. Full documentation at {2}").format(
            phydmslib.__acknowledgments__, phydmslib.__version__,
            phydmslib.__url__),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help='Output file prefix.', type=str)
    parser.add_argument('alignment', help='Existing FASTA file with aligned '
            'codon sequences.', type=ExistingFile)
    parser.add_argument('prefsfiles', help='Existing files with site-specific '
            'amino-acid preferences.', type=ExistingFile, nargs='+')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--raxml', help="Path to RAxML (e.g., 'raxml')")
    group.add_argument('--tree', type=ExistingFile,
             help="Existing Newick file giving input tree.")
    parser.add_argument('--ncpus', default=-1, help='Use this many CPUs; -1 '
            'means all available.', type=int)
    parser.add_argument('--brlen', choices=['scale', 'optimize'],
            default='optimize', help=("How to handle branch lengths: "
            "scale by single parameter or optimize each one"))
    parser.set_defaults(omegabysite=False)
    parser.add_argument('--omegabysite', dest='omegabysite',
            action='store_true', help="Fit omega (dN/dS) for each site.")
    parser.set_defaults(diffprefsbysite=False)
    parser.add_argument('--diffprefsbysite', dest='diffprefsbysite',
            action='store_true', help="Fit differential preferences for "
            "each site.")
    parser.set_defaults(gammaomega=False)
    parser.add_argument('--gammaomega', dest='gammaomega', action=\
            'store_true', help="Fit ExpCM with gamma distributed omega.")
    parser.set_defaults(gammabeta=False)
    parser.add_argument('--gammabeta', dest='gammabeta', action=\
            'store_true', help="Fit ExpCM with gamma distributed beta.")
    parser.set_defaults(noavgprefs=False)
    parser.add_argument('--no-avgprefs', dest='noavgprefs', action='store_true',
            help="No fitting of models with preferences averaged across sites "
            "for ExpCM.")
    parser.set_defaults(randprefs=False)
    parser.add_argument('--randprefs', dest='randprefs', action='store_true',
            help="Include ExpCM models with randomized preferences.")
    parser.add_argument('-v', '--version', action='version', version=
            '%(prog)s {version}'.format(version=phydmslib.__version__))
    return parser


def PhyDMSParser():
    """Returns *argparse.ArgumentParser* for ``phydms`` script."""
    parser = ArgumentParserNoArgHelp(description=('Phylogenetic analysis '
            'informed by deep mutational scanning data. {0} Version {1}. Full'
            ' documentation at {2}').format(phydmslib.__acknowledgments__,
            phydmslib.__version__, phydmslib.__url__),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('alignment', type=ExistingFile,
            help='Existing FASTA file with aligned codon sequences.')
    parser.add_argument('tree', type=ExistingFile,
            help="Existing Newick file giving input tree.")
    parser.add_argument('model', type=ModelOption,
            help=("Substitution model: ExpCM_<prefsfile> or YNGKP_<m> ("
            "where <m> is {0}). For ExpCM, <prefsfile> has first "
            "column labeled 'site' and others labeled by 1-letter "
            "amino-acid code.").format(', '.join(yngkp_modelvariants)))
    parser.add_argument('outprefix', help='Output file prefix.', type=str)
    parser.add_argument('--brlen', choices=['scale', 'optimize'],
            default='optimize', help=("How to handle branch lengths: "
            "scale by single parameter or optimize each one"))
    parser.set_defaults(gammaomega=False)
    parser.add_argument('--gammaomega', action='store_true',
            dest='gammaomega', help="Omega for ExpCM from gamma "
            "distribution rather than single value. To achieve "
            "same for YNGKP, use 'model' of YNGKP_M5.")
    parser.set_defaults(gammabeta=False)
    parser.add_argument('--gammabeta', action='store_true',
            dest='gammabeta', help="Beta for ExpCM from gamma "
            "distribution rather than single value.")
    parser.set_defaults(omegabysite=False)
    parser.add_argument('--omegabysite', dest='omegabysite',
            action='store_true', help="Fit omega (dN/dS) for each site.")
    parser.set_defaults(omegabysite_fixsyn=False)
    parser.add_argument('--omegabysite_fixsyn', dest='omegabysite_fixsyn',
            action='store_true', help="For '--omegabysite', assign all "
            "sites same dS rather than fit for each site.")
    parser.set_defaults(diffprefsbysite=False)
    parser.add_argument('--diffprefsbysite', dest='diffprefsbysite',
            action='store_true', help="Fit differential preferences "
            "for each site.")
    parser.add_argument('--diffprefsprior', default='invquadratic,150,0.5',
            type=diffPrefsPrior, help="Regularizing prior for "
            "'--diffprefsbysite': 'invquadratic,C1,C2' is prior in "
            "Bloom, Biology Direct, 12:1.")
    parser.set_defaults(fitphi=False)
    parser.add_argument('--fitphi', action='store_true', dest='fitphi',
            help='Fit ExpCM phi rather than setting so stationary '
            'state matches alignment frequencies.')
    parser.set_defaults(randprefs=False)
    parser.add_argument('--randprefs', dest='randprefs', action='store_true',
            help="Randomize preferences among sites for ExpCM.")
    parser.set_defaults(avgprefs=False)
    parser.add_argument('--avgprefs', dest='avgprefs', action='store_true',
            help="Average preferences across sites for ExpCM.")
    parser.add_argument('--divpressure', type=ExistingFileOrNone,
            help=("Known diversifying pressure at sites: file with column 1 "
            "= position, column 2 = diversification pressure; columns space-, "
            "tab-, or comma-delimited."))
    parser.add_argument('--ncpus', default=1, type=int,
            help='Use this many CPUs; -1 means all available.')
    parser.add_argument('--fitprefsmethod', choices=[1, 2], default=2,
            help='Implementation to we use when fitting prefs.', type=int)
    parser.add_argument('--ncats', default=4, type=IntGreaterThanOne,
            help='Number of categories for gamma-distribution.')
    parser.add_argument('--minbrlen', type=FloatGreaterThanZero,
            default=phydmslib.constants.ALMOST_ZERO,
            help="Adjust all branch lengths in starting 'tree' to >= this.")
    parser.add_argument('--minpref', default=0.002, type=FloatGreaterThanZero,
            help="Adjust all preferences in ExpCM 'prefsfile' to >= this.")
    parser.add_argument('--seed', type=int, default=1, help="Random number seed.")
    parser.add_argument('--initparams', type=ExistingFile, help="Initialize "
            "model params from this file, which should be format of "
            "'*_modelparams.txt' file created by 'phydms' with this model.")
    parser.set_defaults(profile=False)
    parser.add_argument('--profile', dest='profile', action='store_true',
            help="Profile likelihood maximization, write pstats files. "
            "For code-development purposes.")
    parser.set_defaults(opt_details=False)
    parser.add_argument('--opt_details', dest='opt_details',
            action='store_true', help='Print details about optimization')
    parser.set_defaults(nograd=False)
    parser.add_argument('--nograd', dest='nograd', action='store_true',
            help="Do not use gradients for likelihood maximization.")
    parser.add_argument('-v', '--version', action='version', version=(
            ('%(prog)s {version}'.format(version=phydmslib.__version__))))
    return parser


def PhyDMSTestdivpressureParser():
    """Returns *argparse.ArgumentParser* for ``phdyms_testdivpressure``
        script.
    """
    parser = ArgumentParserNoArgHelp(description="Test different models of "
            "diversifying pressure. {0} Version {1}. Full documentation at {2}"\
            .format(phydmslib.__acknowledgments__, phydmslib.__version__,\
            phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help='Prefix for output files.', type=str)
    parser.add_argument('alignment', help='Existing FASTA file with aligned '
            'codon sequences.', type=ExistingFile)
    parser.add_argument('prefsfile', help='Existing file with site-specific'
            ' amino-acid preferences.', type=ExistingFile)
    parser.add_argument('divpressure', help='List of existing files with '
            'diversifying pressure at each site', type=ExistingFile, nargs='+')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--tree', default=False, type=ExistingFile,
             help="Existing Newick file giving input tree.")
    group.add_argument('--raxml', help="Path to RAxML (http://sco.h-its.org/"
            "exelixis/software.html).", default='raxml')
    parser.add_argument('--randomizations', help = 'Number diversifying '
            'pressure randomizations.', default=0, type=NonNegativeInt)
    parser.add_argument('--ncpus', default=-1, help='Use this many CPUs; '
            '-1 means all available.', type=int)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '
            '{version}'.format(version=phydmslib.__version__))
    return parser


if __name__ == '__main__':
    import doctest
    doctest.testmod()
