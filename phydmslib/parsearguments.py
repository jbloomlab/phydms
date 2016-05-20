"""
Module for parsing arguments.
"""


import sys
import os
import re
import argparse
import phydmslib


# allowed variants of YNGKP models
yngkp_modelvariants = ['M0', 'M1', 'M2', 'M3', 'M7', 'M8']

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

def IntGreaterThanZero(n):
    """If *n* is an integter > 0, returns it, otherwise an error."""
    try:
        n = int(n)
    except:
        raise argparse.ArgumentTypeError("%s is not an integer" % n)
    if n <= 0:
        raise argparser.ArgumentTypeError("%d is not > 0" % n)
    else:
        return n

def FloatGreaterThanEqualToZero(x):
    """If *x* is a float >= 0, returns it, otherwise raises and error.

    >>> print('%.1f' % FloatGreaterThanEqualToZero('1.5'))
    1.5

    >>> print('%.1f' % FloatGreaterThanEqualToZero('-1.1'))
    Traceback (most recent call last):
       ...
    ArgumentTypeError: -1.1 not float greater than or equal to zero
    """
    try:
        x = float(x)
    except:
        raise argparse.ArgumentTypeError("%r not float greater than or equal to zero" % x)
    if x >= 0:
        return x
    else:
        raise argparse.ArgumentTypeError("%r not float greater than or equal to zero" % x)


def FloatGreaterThanOne(x):
    """If *x* is a string for a float > 1, returns it, otherwise an error."""
    x = float(x)
    if x > 1:
        return x
    else:
        raise argparse.ArgumentTypeError("%r not a float greater than one" % x)


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
    """If *fname* is name of an existing file return it, otherwise an error.
    
    *fname* can also be the string 'None', in which case we return *None*."""
    if os.path.isfile(fname):
        return fname
    elif fname.lower() == 'none':
        return None
    else:
        raise argparse.ArgumentTypeError("%s must specify a valid file name or 'None'" % fname)


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


def YNGKPList(modellist):
    """Returns list of YNGKP model variants if *modellist* is valid.

    Removes *M0* if present.
        
    >>> YNGKPList("M0,M3,M7")
    ['M3', 'M7']
    """
    models = [m for m in modellist.split(',') if m and m != 'M0'] # don't count M0 as always included
    if not all([m in yngkp_modelvariants for m in models]):
        raise argparse.ArgumentTypeError("YNGKP model list has invalid entries: {0}".format(modellist))
    if len(models) != len(set(models)):
        raise argparse.ArgumentTypeError("YNGKP model list has duplicated entries: {0}".format(modellist))
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
            raise argparse.ArgumentTypeError("ExpCM_ must be followed by the name of an existing file. You specified the following, which is not an existing file: %s" % fname)
    else:
        raise argparse.ArgumentTypeError("Invalid model")


def PhyDMSAnalyzeSelectionParser():
    """Returns an *argparse.ArgumentParser* for ``phydms_analyzeselection``."""
    parser = ArgumentParserNoArgHelp(description="Visualizes distributions of per-site selection inferred with 'phydms'. Can also highlight information for a subset of sites. %s Version %s. Full documentation at %s" % (phydmslib.__acknowledgments__, phydmslib.__version__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('outprefix', help="Prefix for output files.")
    parser.add_argument('selectionfiles', nargs='+', help="Per-site selection file(s) created by 'phydms' (i.e. '*_omegabysite.txt', '*_stringencybysite.txt', '*_diffprefsbysite.txt').", type=ExistingFile)
    parser.add_argument('--names', nargs='+', help="Name(s) describing type of selection to go with each 'selectionfile'.")
    parser.add_argument('--selectedsites', nargs='+', type=ExistingFile, help="File(s) listing selected sites to display as points. If you give one file, it applies to all 'selectionfiles'. Otherwise list a different file (or 'None') for each file in 'selectionfiles'. Column 1 is site number, which can be followed by # giving notes. Lines beginning with '#' are ignored.")
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

    return parser


def PhyDMSRenumberParser():
    """Returns an *argparse.ArgumentParser* for ``phydms_renumber``."""
    parser = ArgumentParserNoArgHelp(description="Renumber by-site output files from 'phydms'. %s Version %s. Full documentation at %s" % (phydmslib.__acknowledgments__, phydmslib.__version__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('renumberfile', type=ExistingFile, help="Column 1 lists current number for a site; column 2 gives new number. Put 'None' in column 2 if you want a site excluded from renumbered output. Lines beginning with '#' are ignored.")
    parser.add_argument('outprefixes', help="Output prefixes for which we renumber files with appropriate suffixes. Use '*' as a wildcard character.", nargs='+')
    parser.add_argument('--renumberedprefix', default='renumbered', help='Add this prefix followed by underscore to names of renumbered files.')
    return parser



def PhyDMSPlotSelectionParser():
    """Returns an *argparse.ArgumentParser* for ``phydms_plotselection``."""
    parser = ArgumentParserNoArgHelp(description="Visualization of site-specific selection inferred using 'phydms' with an 'ExpCM' model. %s Version %s. Full documentation at %s" % (phydmslib.__acknowledgments__, phydmslib.__version__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('diffprefsbysite', help="A '*_diffprefsbysite.txt' file created by 'phydms' specifying the differential preferences for each site.", type=ExistingFile)
    parser.add_argument('plotfile', help='Name of created PDF file.')
    parser.add_argument('--omegabysite', help="To overlay site-specific omega values, specify a '*_omegabysite.txt' file created by 'phydms'. Sites must match those in 'diffprefs'.", type=ExistingFile)
    parser.add_argument('--stringencybysite', help="To overlay site-specific stringency (beta) values, specify a '*_stringencybysite.txt' file created by 'phydms'. Sites must match those in 'diffprefs'.", type=ExistingFile)
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
    parser.add_argument('--dateseqs', default='None', type=ExistingFile, help="See 'phydms' option of the same name.")
    parser.set_defaults(gammarates=False)
    parser.add_argument('--gammarates', dest='gammarates', action='store_true', help="See 'phydms' option of the same name.")
    parser.add_argument('--avgrandcontrol', default='None', type=ExistingFile, help="Fit average and random controls only for ExpCM using this preference file. Overrides '--no-avgprefs' and '--randprefs'.")
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
    parser = ArgumentParserNoArgHelp(description='Phylogenetic inference using deep mutational scanning data. %s Version %s. Full documentation at %s' % (phydmslib.__acknowledgments__, phydmslib.__version__, phydmslib.__url__), formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('alignment', help='Existing FASTA file with aligned codon sequences.', type=ExistingFile)
    parser.add_argument('tree', help="Existing Newick tree file or 'random' or 'nj'.", type=TreeFile)
    parser.add_argument('model', help="Codon substitution model: YNGKP_<m> where <m> is {0}; or ExpCM_<prefsfile>".format(', '.join(yngkp_modelvariants)), type=ModelOption)
    parser.add_argument('outprefix', help='Prefix for output files.', type=str)
    parser.set_defaults(omegabysite=False)
    parser.add_argument('--omegabysite', dest='omegabysite', action='store_true', help="Fit a different omega (dN/dS) for each site, similar to FEL.")
    parser.set_defaults(stringencybysite=False)
    parser.add_argument('--stringencybysite', dest='stringencybysite', action='store_true', help="Fit a different stringency parameter for each site, only for ExpCM.")
    parser.set_defaults(diffprefsbysite=False)
    parser.add_argument('--diffprefsbysite', dest='diffprefsbysite', action='store_true', help="Infer differential preferences for each site, only for ExpCM.")
    parser.set_defaults(infertopology=False)
    parser.add_argument('--infertopology', dest='infertopology', action='store_true', help="Infer topology starting from 'tree'; otherwise topology fixed to 'tree'. Requires YNGKP_M0 model.")
    parser.set_defaults(gammarates=False)
    parser.add_argument('--gammarates', dest='gammarates', action='store_true', help="Draw the substitution rate from 4-category discrete gamma distribution with optimized shape parameter.")
    parser.add_argument('--ncats', default=6, help="Number of beta-distributed categories for YNGKP M7 and M8.", type=IntGreaterThanZero)
    parser.add_argument('--dateseqs', type=ExistingFile, help="Perform least-squares sequence dating. Specify file with column 1 = date, column 2 = sequence header; columns are space delimited.")
    parser.add_argument('--ncpus', default=1, help='Use this many CPUs; -1 means all available.', type=int)
    parser.set_defaults(omegabysite_fixsyn=False)
    parser.add_argument('--omegabysite_fixsyn', dest='omegabysite_fixsyn', action='store_true', help="For '--omegabysite', assign all sites same synonymous rate rather than fitting a different one for each site.")
    # comment out this option as the 'invquad' prior makes vastly more sense than 'dirichlet'
    parser.set_defaults(diffprefsprior='invquad')
    #parser.add_argument('--diffprefsprior', default='invquad', choices=['invquad', 'dirichlet'], help="Prior on diff preferences for '--diffprefsbysite'.")
    parser.set_defaults(randprefs=False)
    parser.add_argument('--randprefs', dest='randprefs', action='store_true', help="Randomize preferences among sites for ExpCM.")
    parser.set_defaults(avgprefs=False)
    parser.add_argument('--avgprefs', dest='avgprefs', action='store_true', help="Average preferences across sites for ExpCM.")
    parser.set_defaults(fixbrlen=False)
    parser.add_argument('--fixbrlen', dest='fixbrlen', action='store_true', help="Fix branch lengths to those of initial 'tree'. Consider using '--addrateparameter' too.")
    parser.add_argument('--diffprefconc', help="Parameters determining the concentration of the regularizing prior over the differential preferences for '--diffprefsbysite'. Larger values favor smaller diff prefs.", type=FloatGreaterThanEqualToZero, nargs=2, metavar=('C1', 'C2'), default=[150, 0.5])
    parser.set_defaults(addrateparameter=False)
    parser.add_argument('--addrateparameter', dest='addrateparameter', action='store_true', help="Add parameter scaling substitution rate. Only allowed with '--fixbrlen'.")
    parser.set_defaults(fitF3X4=False)
    parser.add_argument('--fitF3X4', dest='fitF3X4', action='store_true', help='Fit F3X4 frequencies for YNGKP; otherwise use empirical')
    parser.add_argument('--minbrlen', default=1e-6, help='Min branch length for starting tree.', type=FloatGreaterThanZero)
    parser.add_argument('--seed', type=int, default=1, help="Random number seed.")
    parser.set_defaults(no_optimize=False)
    parser.add_argument('--no_optimize', dest='no_optimize', action='store_true', help="Don't optimize tree or model; use values in existing files from previous run with same 'outprefix'.")
    parser.set_defaults(useLog=False)
    parser.add_argument('--useLog', dest='useLog', action='store_true', help="Use logarithms in likelihood calculations.")
    parser.add_argument('--debugsite', type=int, help="For debugging: fit site-specific selection to only this site.")
    parser.set_defaults(recursion='S')
    # comment out this option for now as the double recursion seems not to work properly
    #parser.add_argument('--recursion', choices=['S', 'D'], default='S', help='Likelihood recursion for Bio++.')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=phydmslib.__version__))
    return parser


if __name__ == '__main__':
    import doctest
    doctest.testmod()
