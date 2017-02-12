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


def ExistingFile(fname):
    """If *fname* is name of an existing file return it, otherwise an error.
    
    *fname* can also be the string 'None', in which case we return *None*."""
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


def YNGKPList(modellist):
    """Returns list of YNGKP model variants if *modellist* is valid.

    Removes *M0* if present.
        
    >>> YNGKPList("M0,M5")
    ['M5']
    """
    models = [m for m in modellist.split(',') if m and m != 'M0'] # don't count M0 as always included
    if not all([m in yngkp_modelvariants for m in models]):
        raise ValueError("YNGKP model list has invalid entries: {0}".format(modellist))
    if len(models) != len(set(models)):
        raise ValueError("YNGKP model list has duplicated entries: {0}".format(modellist))
    return models


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


def PhyDMSAnalyzeSelectionParser():
    """Returns an *argparse.ArgumentParser* for ``phydms_analyzeselection``."""
    parser = ArgumentParserNoArgHelp(description="Visualizes distributions of per-site selection inferred with 'phydms'. Can also highlight information for a subset of sites. %s Version %s. Full documentation at %s" % (phydmslib.__acknowledgments__, phydmslib.__version__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help="Prefix for output files.")
    parser.add_argument('selectionfiles', nargs='+', help="Per-site selection file(s) created by 'phydms' (i.e. '*_omegabysite.txt', '*_stringencybysite.txt', '*_diffprefsbysite.txt').", type=ExistingFile)
    parser.add_argument('--names', nargs='+', help="Name(s) describing type of selection to go with each 'selectionfile'.")
    parser.add_argument('--selectedsites', nargs='+', type=ExistingFileOrNone, help="File(s) listing selected sites to display as points. If you give one file, it applies to all 'selectionfiles'. Otherwise list a different file (or 'None') for each file in 'selectionfiles'. Column 1 is site number, which can be followed by # giving notes. Lines beginning with '#' are ignored.")
    parser.set_defaults(labelselectedsites=False)
    parser.add_argument('--labelselectedsites', action='store_true', dest='labelselectedsites', help="Do we use a unique labeled point on violin plots for each site in '--selectedsites'?")
    parser.add_argument('--fdr', type=float, default=0.05, help="False discovery rate for 'omega' and 'stringency'. Benjamini-Hochberg FDR computed separately for values > and < 1.")
    parser.add_argument('--maxlog10p', type=FloatGreaterThanZero, default=5, help="For 'omega' and 'stringency' violin plots, if log10 P-value has magnitude > this, instead plot as this. Also is y-limits for these plots.")
    parser.set_defaults(groupbyname=False)
    parser.add_argument('--groupbyname', action='store_true', dest='groupbyname', help="Group selections with same first word in name specified by '--names'. Grouped selections must have same sites specified by '--selectedsites'.")
    parser.add_argument('--diffprefsline', help="Draw cutoff line on diffprefs plot: 'lowestpeak' indicates the lowest peak value for any selection file, otherwise specify a number between 0 and 1.")
    parser.set_defaults(nolegend=False)
    parser.add_argument('--nolegend', action='store_true', dest='nolegend', help="Don't place a legend even on plates with labeled selected sites.")
    parser.set_defaults(dNdSlabel=False)
    parser.add_argument('--dNdSlabel', action='store_true', dest='dNdSlabel', help="Label omega-by-site plots with 'dN/dS' rather than 'omega_r'.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=phydmslib.__version__))
    return parser


def PhyDMSRenumberParser():
    """Returns an *argparse.ArgumentParser* for ``phydms_renumber``."""
    parser = ArgumentParserNoArgHelp(description="Renumber by-site output files from 'phydms'. %s Version %s. Full documentation at %s" % (phydmslib.__acknowledgments__, phydmslib.__version__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('renumberfile', type=ExistingFile, help="Column 1 lists current number for a site; column 2 gives new number. Put 'None' in column 2 if you want a site excluded from renumbered output. Lines beginning with '#' are ignored.")
    parser.add_argument('outprefixes', help="Output prefixes for which we renumber files with appropriate suffixes. Use '*' as a wildcard character.", nargs='+')
    parser.add_argument('--renumberedprefix', default='renumbered', help='Add this prefix followed by underscore to names of renumbered files.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=phydmslib.__version__))
    return parser


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


def PhyDMSPlotSelectionParser():
    """Returns an *argparse.ArgumentParser* for ``phydms_plotselection``."""
    parser = ArgumentParserNoArgHelp(description="Visualization of site-specific selection inferred using 'phydms' with an 'ExpCM' model. %s Version %s. Full documentation at %s" % (phydmslib.__acknowledgments__, phydmslib.__version__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('diffprefsbysite', help="A '*_diffprefsbysite.txt' file created by 'phydms' specifying the differential preferences for each site.", type=ExistingFile)
    parser.add_argument('plotfile', help='Name of created PDF file.')
    parser.add_argument('--omegabysite', help="To overlay site-specific omega values, specify a '*_omegabysite.txt' file created by 'phydms'. Sites must match those in 'diffprefs'.", type=ExistingFileOrNone)
    parser.add_argument('--stringencybysite', help="To overlay site-specific stringency (beta) values, specify a '*_stringencybysite.txt' file created by 'phydms'. Sites must match those in 'diffprefs'.", type=ExistingFileOrNone)
    parser.add_argument('--nperline', type=IntGreaterThanZero, default=70, help="Number of sites per line in plot.")
    parser.add_argument('--numberevery', type=IntGreaterThanZero, default=10, help="Number sites at this interval.")
    parser.add_argument('--diffprefheight', type=FloatGreaterThanZero, default=1.0, help="Height of differential preferences logo stacks in each direction. If using '--updiffprefheight' then the height may be higher than this.")
    parser.set_defaults(updiffprefheight=False)
    parser.add_argument('--updiffprefheight', dest='updiffprefheight', action='store_true', help="Automatically increase '--diffprefheight' to make it exceed max differential preferences.")
    parser.add_argument('--minP', type=FloatGreaterThanZero, default=1e-4, help="Minimum plotted P-value for omega and stringency by site.")
    parser.add_argument('--colormap', type=str, default='jet', help='Colormap for amino-acid hydrophobicity. Must specify a valid ``pylab`` color map')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=phydmslib.__version__))
    return parser


def PhyDMSComprehensiveParser():
    """Returns *argparse.ArgumentParser* for ``phdyms_comprehensive`` script."""
    parser = ArgumentParserNoArgHelp(description="Comprehensive phylogenetic model comparison and detection of selection using deep mutational scanning data. This program runs 'phydms' to infer a tree topology, then compare substitution models, then detect selection at each site. %s Version %s. Full documentation at %s" % (phydmslib.__acknowledgments__, phydmslib.__version__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help='Prefix for output files.', type=str)
    parser.add_argument('alignment', help='Existing FASTA file with aligned codon sequences.', type=ExistingFile)
    parser.add_argument('prefsfiles', help='Existing files with site-specific amino-acid preferences.', type=ExistingFile, nargs='+')
    parser.add_argument('--treetopology', help='Fix tree to this Newick topology.', default=None)
    parser.add_argument('--ncpus', default=-1, help='Use this many CPUs; -1 means all available.', type=int)
    parser.set_defaults(noomegabysite=False)
    parser.add_argument('--no-omegabysite', dest='noomegabysite', action='store_true', help="No fitting of site-specific omegas.")
    parser.set_defaults(nostringencybysite=False)
    parser.add_argument('--no-stringencybysite', dest='nostringencybysite', action='store_true', help="No fitting of site-specific stringency for ExpCM.")
    parser.set_defaults(nodiffprefsbysite=False)
    parser.add_argument('--no-diffprefsbysite', dest='nodiffprefsbysite', action='store_true', help="No fitting of differential preferences for ExpCM.")
    parser.set_defaults(noavgprefs=False)
    parser.add_argument('--no-avgprefs', dest='noavgprefs', action='store_true', help="No fitting of models with preferences averaged across sites for ExpCM.")
    parser.add_argument('--dateseqs', default='None', type=ExistingFileOrNone, help="See 'phydms' option of the same name.")
    parser.set_defaults(gammarates=False)
    parser.add_argument('--gammarates', dest='gammarates', action='store_true', help="See 'phydms' option of the same name.")
    parser.add_argument('--avgrandcontrol', default='None', type=ExistingFileOrNone, help="Fit average and random controls only for ExpCM using this preference file. Overrides '--no-avgprefs' and '--randprefs'.")
    parser.set_defaults(omegabysite_fixsyn=False)
    parser.add_argument('--omegabysite_fixsyn', action='store_true', dest='omegabysite_fixsyn', help="See 'phydms' option of the same name.")
    parser.set_defaults(useLog=False)
    parser.add_argument('--useLog', action='store_true', dest='useLog', help="See 'phydms' option of the same name.")
    parser.add_argument('--yngkp', default='M8', help="YNGKP models to use in addition to M0. Should be comma-separated list of models from the following: {0}".format(','.join([m for m in yngkp_modelvariants if m != 'M0'])), type=YNGKPList)
    parser.add_argument('--ncats', default=6, help="Number of beta-distributed categories for YNGKP M7 and M8.", type=IntGreaterThanZero)
    parser.add_argument('--diffprefconc', help="Parameters determining the concentration of the regularizing prior over the differential preferences for '--diffprefsbysite'. Defaults to value of the same param for 'phydms'.", type=FloatGreaterThanEqualToZero, nargs=2, metavar=('C1', 'C2'))
    parser.set_defaults(randprefs=False)
    parser.add_argument('--randprefs', dest='randprefs', action='store_true', help="Include ExpCM models with randomized preferences.")
    parser.set_defaults(fitF3X4=False)
    parser.add_argument('--fitF3X4', dest='fitF3X4', action='store_true', help='Fit F3X4 frequencies for YNGKP; otherwise use empirical')
    parser.set_defaults(use_existing=False)
    parser.add_argument('--use_existing', action='store_true', dest='use_existing', help="Use existing 'phydms' output for a model if it exists. BE CAREFUL: no checks are performed to ensure calling options the same.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=phydmslib.__version__))
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
    parser.add_argument('--brlen', choices=['scale', 'optimize', 'fix'],
            default='scale', help=("How to handle branch lengths: scale "
            "by single parameter; optimize each one; fix to values "
            "in 'tree'.")) 
    parser.set_defaults(gammaomega=False)
    parser.add_argument('--gammaomega', action='store_true', 
            dest='gammaomega', help="Omega for ExpCM from gamma "
            "distribution rather than single value. To achieve "
            "same for YNGKP, use 'model' of YNGKP_M5.")
    parser.set_defaults(fitphi=False)
    parser.add_argument('--fitphi', action='store_true', dest='fitphi',
            help='Fit ExpCM phi rather than setting so stationary '
            'state matches alignment frequencies.')
    parser.set_defaults(omegabysite=False)
    parser.add_argument('--omegabysite', dest='omegabysite', 
            action='store_true', help="Fit omega (dN/dS) for each site.")
    parser.set_defaults(omegabysite_fixsyn=False)
    parser.add_argument('--omegabysite_fixsyn', dest='omegabysite_fixsyn', 
            action='store_true', help="For '--omegabysite', assign all "
            "sites same dS rather than for each site.")
    parser.set_defaults(randprefs=False)
    parser.add_argument('--randprefs', dest='randprefs', action='store_true', 
            help="Randomize preferences among sites for ExpCM.")
    parser.set_defaults(avgprefs=False)
    parser.add_argument('--avgprefs', dest='avgprefs', action='store_true', 
            help="Average preferences across sites for ExpCM.")
    parser.add_argument('--divpressure', type=ExistingFileOrNone, 
            help=("Known diversifying pressure at sites: file with column 1 "
            "= position, column 2 = diversification pressure; columns space "
            "delimited."))
    parser.add_argument('--ncpus', default=1, type=int,
            help='Use this many CPUs; -1 means all available.')
    parser.add_argument('--ncats', default=4, type=IntGreaterThanOne,
            help='Number of categories for gamma-distributed omega.')
    parser.add_argument('--minbrlen', type=FloatGreaterThanZero,
            default=phydmslib.constants.ALMOST_ZERO, 
            help="Adjust all branch lengths in starting 'tree' to >= this.")
    parser.add_argument('--minpref', default=0.005, type=FloatGreaterThanZero,
            help="Adjust all preferences in ExpCM 'prefsfile' to >= this.")
    parser.add_argument('--seed', type=int, default=1, help="Random number seed.")
    parser.set_defaults(profile=False)
    parser.add_argument('--profile', dest='profile', action='store_true', 
            help="Profile likelihood maximization, write pstats files. "
            "For code-development purposes.")
    parser.add_argument('-v', '--version', action='version', version=(
            ('%(prog)s {version}'.format(version=phydmslib.__version__))))
    return parser
    
    
def PhyDMSTestdivpressureParser():
    """Returns *argparse.ArgumentParser* for ``phdyms_testdivpressure`` script."""
    parser = ArgumentParserNoArgHelp(description="Test different models of diversifying pressure. This program runs 'phydms' and compares model log likelihoods with and without diversifying pressures. %s Version %s. Full documentation at %s" % (phydmslib.__acknowledgments__, phydmslib.__version__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help='Prefix for output files.', type=str)
    parser.add_argument('alignment', help='Existing FASTA file with aligned codon sequences.', type=ExistingFile)
    parser.add_argument('prefsfile', help='Existing file with site-specific amino-acid preferences.', type=ExistingFile)
    parser.add_argument('divpressure', help='List of existing files with diversifying pressure at each site', type=ExistingFile, nargs='+')
    parser.add_argument('--treetopology', help='Fix tree to this Newick topology.', default=None)
    parser.add_argument('--randomizations', help = 'Number diversifying pressure randomizations.', default=0, type=NonNegativeInt)
    parser.add_argument('--ncpus', default=-1, help='Use this many CPUs; -1 means all available.', type=int)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=phydmslib.__version__))
    return parser


if __name__ == '__main__':
    import doctest
    doctest.testmod()
