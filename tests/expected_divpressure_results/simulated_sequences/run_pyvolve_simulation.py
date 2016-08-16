"""Simulates alignment using Bloom2016 model with ``pyvolve``.

Run script to get help message with info on command::

    python run_simulation.py

Written by Jesse Bloom

Edited by Hugh Haddox, May-29-2016

Edited by Sarah Hilton, August-15-2016
"""


import re
import os
import sys
import math
import random
import argparse
import Bio.SeqIO
import pyvolve
import dms_tools.file_io


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""
    def error(self, message):
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def ParseArguments():
    """Argument parser for script."""
    parser = ArgumentParserNoArgHelp(
            description='Simulates alignment according to Bloom2016 model with pyvolve.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            )
    parser.add_argument('prefs', help='File containing preferences.')
    parser.add_argument('tree', help='File containing tree.')
    parser.add_argument('modelparams', help='File containing model params from phydms.')
    parser.add_argument('simulatedalignment', help='Name of created alignment file')
    parser.add_argument('--seed', type=int, default=1, help='Random number seed')
    return parser


def ReadParams(paramsfile):
    """Reads parameters of Bloom2016 model from phydms modelparams file.
    
    Returned dictionary has the following keys:
    *phiA*, *phiC*, *phiG*, *phiT*, *kappa*, *stringencyparameter*
    """
    possibleParams = ['stringencyparameter', 'kappa', 'scalerate','diversifyingomegaA','diversifyingsitesA','diversifyingomegaB','diversifyingsitesB'] +\
                         ['phi{0}'.format(nt) for nt in ['A', 'T', 'G', 'C']]
    params = {}
    for param in possibleParams:
        params[param] = None
    with open(paramsfile) as f:
        paramstext = f.read()
    for param in params.keys():
        key = param.split('.')[-1]
        if not key.startswith('diversifyingsites'):
        	if bool(re.search('{0} = (\d+\.\d+)'.format(param), paramstext)):
        	    params[key] = float(re.search('{0} = (\d+\.\d+)'.format(param), paramstext).group(1))
        else:
            if bool(re.search('{0} = (\d+.+)'.format(param), paramstext)):
        	    params[key] = [int(x) for x in re.search('{0} = (\d+.+)'.format(param), paramstext).group(1).split(',')]
    print params
    return params

def main():
    """Main body of script."""
    codons = pyvolve.genetics.Genetics().codons
    codon_dict = pyvolve.genetics.Genetics().codon_dict
    pyrims = pyvolve.genetics.Genetics().pyrims
    purines = pyvolve.genetics.Genetics().purines

    args = vars(ParseArguments().parse_args())
    print("Read the following command line arguments:")
    print("\n\t{0}".format("\n\t".join(["{0} = {1}".format(key, value) for (key, value) in args.items()])))

    print("\nPerforming simulation with pyvolve version {0}".format(pyvolve.__version__))

    print("\nReading model params from {0}".format(args['modelparams']))
    params = ReadParams(args['modelparams'])
    for (param, paramvalue) in params.items():
        print("The value of {0} is {1}".format(param, paramvalue))
        
    print("\nReading preferences from {0}".format(args['prefs']))
    tup = dms_tools.file_io.ReadPreferences(args['prefs'])
    (sites, pis) = (tup[0], tup[2])
    print("\nRead amino-acid preferences for {0} sites".format(len(pis)))

    tree = pyvolve.read_tree(file=args['tree'])

    # create models for simulation
    partitions = []
    for r in sites:
        if params['diversifyingsitesA'] and (int(r) in params['diversifyingsitesA']):
            omega = params['diversifyingomegaA']
            print r,omega
        elif params['diversifyingsitesB'] and (int(r) in params['diversifyingsitesB']):
            omega = params['diversifyingomegaB']
            print r,omega
        else:
            omega = 1.0
        matrix = [] # matrix[x][y] is rate of substitution from x to y
        for (xi, x) in enumerate(codons):
            row = []
            for (yi, y) in enumerate(codons):
                ntdiffs = [(x[j], y[j]) for j in range(3) if x[j] != y[j]]
                if len(ntdiffs) == 0:
                    assert x == y
                    row.append(0) # will later be adjusted to make row sum to zero
                elif len(ntdiffs) > 1:
                    # multi-nucleotide codon change
                    row.append(0)
                else:
                    # single nucleotide change
                    (xnt, ynt) = ntdiffs[0]
                    if (xnt in purines) == (ynt in purines):
                        # transition
                        qxy = params['kappa'] * params['phi{0}'.format(ynt)]
                    else:
                        # transversion
                        qxy = params['phi{0}'.format(ynt)]
                    (xaa, yaa) = (codon_dict[x], codon_dict[y])
                    if xaa == yaa:
                        fxy = 1.0
                    else:
                        pix = pis[r][xaa]**params['stringencyparameter']
                        piy = pis[r][yaa]**params['stringencyparameter']
                        if abs(pix - piy) < 1e-6:
                            fxy = omega
                        else:
                            fxy = omega * math.log(piy / pix) / (1.0 - pix / piy)
                    row.append(qxy * fxy * params['scalerate'])
            assert len(row) == len(codons)
            row[xi] = -sum(row)
            matrix.append(row)          
        model = pyvolve.Model("custom", {"matrix":matrix})
        partitions.append(pyvolve.Partition(models=model, size=1))

    print("\nSimulating evolution, writing to {0}...".format(args['simulatedalignment']))
    basename = os.path.splitext(args['simulatedalignment'])[0]
    evolver = pyvolve.Evolver(partitions=partitions, tree=tree)
    evolver(
            seqfile=args['simulatedalignment'],
            infofile='{0}_infofile.txt'.format(basename),
            ratefile='{0}_ratefile.txt'.format(basename),
            )
    print("Finished simulation")

    uniqueseqs = set([])
    uniquealignment = []
    ninitial = 0
    for seq in Bio.SeqIO.parse(args['simulatedalignment'], 'fasta'):
        ninitial += 1
        seqstr = str(seq.seq)
        if seqstr not in uniqueseqs:
            uniqueseqs.add(seqstr)
            uniquealignment.append(seq)
    print("\nAfter removing redundant sequences, we have shrunk {0} from {1} to {2} sequences".format(args['simulatedalignment'], ninitial, len(uniquealignment)))
    Bio.SeqIO.write(uniquealignment, args['simulatedalignment'], 'fasta')


if __name__ == '__main__':
    main()
