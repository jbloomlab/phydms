"""Module for plotting."""


import math
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt


def PlotSignificantOmega(plotfile, models, ngt, nlt, nsites, fdr, usetex=True):
    """Makes a PDF plot of the number of sites with significant omega.

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
    assert len(models) == len(ngt) == len(nlt)
    plt.rc('font', size=15)
    plt.rc('text', usetex=usetex)
    omegacategories = {'< 1':'bo-', '> 1':'rs-'}
    handles = []
    labels = []
    for (cat, style) in omegacategories.items():
        xs = range(len(models))
        ys = {'< 1':nlt, '> 1':ngt}[cat]
        handle = plt.plot(xs, ys, style, markersize=11)
        handles.append(handle[0])
        if usetex:
            labels.append('$\omega_r %s$' % cat)
        else:
            labels.append('omega %s' % cat)
    plt.xlim(-0.25, len(models) - 0.75)
    plt.xticks(range(len(models)), [modelnames[model] for model in models], fontsize=17)
    plt.locator_params(axis='y', bins=4)
    plt.ylabel('sites out of %d (FDR %.2f)' % (len(sites), fdr_alpha), fontsize=17)
    plt.legend(handles, labels, 'upper right', numpoints=1)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
