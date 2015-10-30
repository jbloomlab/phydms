"""Module for plotting."""


import os
import math
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import dms_tools.utils


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
        xs = [x for x in range(len(models))]
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
    plt.xticks(xs, models, fontsize=32)
    plt.locator_params(axis='y', bins=4)
    plt.ylabel('sites out of %d (FDR %.2f)' % (nsites, fdr), fontsize=30)
    plt.legend(handles, labels, loc='upper right', numpoints=1, fontsize=34, borderaxespad=0)
    plt.savefig(plotfile, bbox_inches='tight')
    plt.clf()
    plt.close()



def SelectionViolinPlot(plotfile, ylabel, models, yvalues, symmetrizey, hlines=None, points=None, pointmarkercolor='or', usetex=True, legend=False, fixymin=None, fixymax=None):
    """Creates violin plot showing distribution of selection and significant sites.

    Calling arguments:

    *plotfile* : name of PDF plot to create.

    *ylabel* : ylabel for the plot.

    *models* : list of models for which we create violin plots.

    *yvalues* : list of the same length as *models*, each entry is a list
    of the Y-values (such as P-values).

    *symmetrizey* : make y-axis symmetric around zero?

    *hlines* : if not *None*, list of the same length as *models* with each entry
    a list giving y-value for where we draw horizontal lines for that *model*.

    *points* : if not *None*, list of the same length as *models* with each entry
    a list giving the y-value for points to be placed for that *model*.

    *pointmarkercolor* : specifies marker and color of points in *points*. 
    Should either a length-two string giving marker and color for all points
    (such as *or* for circles, red) or a list of lists of the same length
    as *points* witch each entry specifying the marker and color for that point.

    *usetex* : use LaTex formatting of strings?

    *legend* : Create a legend with names of points specified by *pointmarkercolor*?
    This option is only valid if *pointmarkercolor* is a list. If it is not *False*,
    it should be the tuple *(markercolors, names)*. In this tuple,
    *markercolors* and *names* are lists of the same
    length, with *markercolors* being a list of marker / color (e.g. *or*
    for circles, red) and *names* being a list of the string corresponding
    to each marker / color.

    *fixymin* : if not *None*, the y-minimum is fixed to this value.

    *fixymax* : if not *None*, the y-maximum is fixed to this value.
    """
    alpha = 0.55 # transparency for points
    markersize = 25 # size of points
    markerlw = 0.6 # line width for makers
    assert os.path.splitext(plotfile)[1].lower() == '.pdf', "plotfile %s does not end with extension '.pdf'"
    assert len(models) == len(yvalues) >= 1
    assert not (legend and not isinstance(pointmarkercolor, list)), "You can only use legend if pointmarkercolor is a list"
    assert not (isinstance(pointmarkercolor, list) and not points), "It doesn't make sense to set pointmarkercolor to a list if you're not using points"
    plt.rc('font', size=12)
    plt.rc('text', usetex=usetex)
    lmargin = 0.7
    tmargin = 0.1
    bmargin = 0.4
    if legend:
        rmargin = 0.95
    else:
        rmargin = 0.1
    (height, widthper) = (2.5, 1.75)
    totwidth = len(models) * widthper + lmargin + rmargin
    totheight = height + tmargin + bmargin
    plt.figure(figsize=(totwidth, totheight))
    plt.axes([lmargin / totwidth, bmargin / totheight, 1.0 - (lmargin + rmargin) / totwidth, 1.0 - (tmargin + bmargin) / totheight])
    plt.ylabel(ylabel, fontsize=16)
    xs = [x for x in range(len(models))]
    violinwidth = 0.75
    plt.violinplot(yvalues, xs, widths=violinwidth, showextrema=False)
    plt.xlim(-1.2 * violinwidth / 2.0, len(models) - 1 + 1.2 * violinwidth / 2.0)
    if hlines:
        assert len(hlines) == len(models)
        line_ys = []
        line_xmins = []
        line_xmaxs = []
        for i in range(len(models)):
            line_ys += hlines[i]
            line_xmins += [i - violinwidth / 2.0] * len(hlines[i])
            line_xmaxs += [i + violinwidth / 2.0] * len(hlines[i])
        plt.hlines(line_ys, line_xmins, line_xmaxs, colors='b')
    if symmetrizey:
        (ymin, ymax) = plt.ylim()
        ymax = 1.05 * max(abs(ymin), abs(ymax))
        ymin = -ymax
    else:
        (ymin, ymax) = plt.ylim()
    if points:
        assert len(points) == len(models)
        if isinstance(pointmarkercolor, str):
            assert len(pointmarkercolor) == 2
            color = pointmarkercolor[1]
            marker = pointmarkercolor[0]
        else:
            color = []
            marker = []
            assert len(pointmarkercolor) == len(points), "len(pointmarkercolor) = %d; len(points) = %d" % (len(pointmarkercolor), len(points))
        point_xs = []
        point_ys = []
        for i in range(len(models)):
            (model_xs, model_ys) = SmartJitter(points[i], yspace=(ymax - ymin) / 25., xspace=0.07, xcenter=i)
            point_xs += model_xs
            point_ys += model_ys
            if not isinstance(pointmarkercolor, str):
                imarkercolor = pointmarkercolor[i]
                assert len(imarkercolor) == len(points[i]), "pointmarkercolor and points have length mismatch for %d" % i
                color += [x[1] for x in imarkercolor]
                marker += [x[0] for x in imarkercolor]
        if isinstance(pointmarkercolor, str):
            plt.scatter(point_xs, point_ys, s=markersize, c=color, marker=marker, alpha=alpha, lw=markerlw)
        else:
            assert len(color) == len(marker) == len(point_xs)
            for (x, y, c, m) in zip(point_xs, point_ys, color, marker):
                plt.scatter(x, y, s=markersize, c=c, marker=m, alpha=alpha, lw=markerlw)
    if fixymin not in [None, False]:
        ymin = fixymin
    if fixymax not in [None, False]:
        ymax = fixymax
    assert ymin < ymax
    plt.ylim(ymin, ymax)
    plt.xticks(xs, models, fontsize=16)
    if legend:
        (markercolors, legendnames) = legend
        assert len(markercolors) == len(legendnames)
        handles = [matplotlib.lines.Line2D([0], [0], marker=marker, color=color, markersize=markersize, alpha=alpha, lw=markerlw, linestyle='None') for (marker, color) in markercolors]
        # put in natural sort order
        assert len(set(legendnames)) == len(legendnames), "Duplicate legendnames entry"
        sortedlegendnames = list(legendnames)
        dms_tools.utils.NaturalSort(sortedlegendnames)
        sortedhandles = [None] * len(sortedlegendnames)
        for (handle, name) in zip(handles, legendnames):
            sortedhandles[sortedlegendnames.index(name)] = handle
        assert None not in sortedhandles
        legend = plt.figlegend(sortedhandles, sortedlegendnames, loc='upper right', fontsize=13, numpoints=1, title='\\bf{sites}', markerscale=0.25, handlelength=0.7, handletextpad=0.25)
        plt.setp(legend.get_title(), fontsize=14)
    plt.savefig(plotfile)
    plt.clf()
    plt.close()


def SmartJitter(ys, yspace, xspace, xcenter):
    """Smartly horizontally spaces points with assigned y-values to decrease overlap.

    Divides y-axis into bins, in each bin spaces points horizontally with the largest
    y-value point in the center and lower y-value points spreading out horizontally.

    *ys* is a list of y-coordinates of the points.

    *yspace* is the spacing between y-axis bins.

    *xspace* is the spacing between points in the same bin on the x-axis.

    *xcenter* is the center for the y-axis.

    The return value is *(smart_xs, smart_ys)*, which is a list with the
    x and y values of the spaced points. 
    """
    assert yspace > 0 and xspace > 0
    (ymin, ymax) = (min(ys), max(ys))
    (ymin, ymax) = (ymin - yspace / 10. - 1e-5, ymax + yspace / 10. + 1e-5)
    assigned = [False] * len(ys) 
    smart_xs = []
    smart_ys = []
    indices = [] # indices[j] is the index of the point in smart_ys[j] in the original xs, ys
    binymin = ymin
    while not all(assigned):
        assert binymin <= ymax
        binymax = binymin + yspace
        yindices = [iy for (iy, y) in enumerate(ys) if binymin <= y < binymax]
        binymin += yspace
        assert all([not assigned[iy] for iy in yindices]), "Assigned point to duplicate bin"
        binys = []
        for iy in yindices:
            assigned[iy] = True
            binys.append((ys[iy], iy))
        if not binys:
            continue
        # make centeredbinys so that largest value is in middle
        binys.sort()
        binys.reverse()
        centeredbinys = [binys[0]]
        before = True
        for y in binys[1 : ]:
            if before:
                centeredbinys.insert(0, y)
            else:
                centeredbinys.append(y)
            before = not before
        assert len(centeredbinys) == len(binys)
        smart_ys += [tup[0] for tup in centeredbinys]
        indices += [tup[1] for tup in centeredbinys]
        xmin = xcenter - xspace * (len(binys) - 1) / 2.0
        smart_xs += [xmin + xspace * i for i in range(len(binys))]
    assert all(assigned), "Failed to assign all points to bins"
    assert set(indices) == set([i for i in range(len(ys))]), "Failed to assign all unique indices"
    assert len(smart_xs) == len(smart_ys) == len(ys) == len(assigned)
    reindexed_smart_xs = [None] * len(smart_xs)
    reindexed_smart_ys = [None] * len(smart_ys)
    for (i, j) in enumerate(indices):
        reindexed_smart_xs[j] = smart_xs[i]
        reindexed_smart_ys[j] = smart_ys[i]
    return (reindexed_smart_xs, reindexed_smart_ys)



if __name__ == '__main__':
    import doctest
    doctest.testmod()
