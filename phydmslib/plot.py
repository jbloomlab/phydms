"""Module for plotting."""


import os
import math
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt


def PlotSignificantOmega(plotfile, models, ngt, nlt, nsites, fdr, usetex=True):
    """Makes a PDF slopegraph of the number of sites with significant omega.

    *plotfile* : name of created PDF.

    *models* : list of model names as strings

    *ngt* : list of same length as *models* giving number of sites with
    *omega* > 1.

    *nlt* : list of same length as *models* giving number of sites with
    *omega* < 1.

    *nsites* : total number of sites (used for y-axis label).

    *fdr* : false discovery rate (used for y-axis label).

    *usetex* : use LaTex formatting of strings?
    """
    assert os.path.splitext(plotfile)[1].lower() == '.pdf', "plotfile %s does not end with extension '.pdf'"
    assert len(models) == len(ngt) == len(nlt)
    plt.rc('font', size=25)
    plt.rc('text', usetex=usetex)
    omegacategories = {'< 1':'bo-', '> 1':'rs-'}
    handles = []
    labels = []
    ymax = 1
    for (cat, style) in omegacategories.items():
        xs = range(len(models))
        ys = {'< 1':nlt, '> 1':ngt}[cat]
        ymax = max(ymax, max(ys))
        handle = plt.plot(xs, ys, style, markersize=22, linewidth=3)
        handles.append(handle[0])
        if usetex:
            labels.append('$\omega_r %s$' % cat)
        else:
            labels.append('omega %s' % cat)
    plt.xlim(-0.25, len(models) - 0.75)
    plt.ylim(0, int(1.02 * ymax + 1))
    plt.xticks(range(len(models)), models, fontsize=32)
    plt.locator_params(axis='y', bins=4)
    plt.ylabel('sites out of %d (FDR %.2f)' % (nsites, fdr), fontsize=30)
    plt.legend(handles, labels, loc='upper right', numpoints=1, fontsize=34, borderaxespad=0)
    plt.savefig(plotfile, bbox_inches='tight')
    plt.clf()
    plt.close()



if __name__ == '__main__':
    import doctest
    doctest.testmod()
