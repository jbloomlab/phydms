"""Module for plotting."""


import os
import math
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import dms_tools.utils


def SplitText(text, maxchars):
    r"""Splits text by word and then breaks long words.

    *text* is the text to split, *maxchars* is the max characters
    per line. First splits every word onto a new line, then 
    splits long words.

    >>> SplitText('lactamase ExpCM', maxchars=5)
    'lact-\namase\nExpCM'
    """
    assert maxchars > 1, 'maxchars must be > 1'
    newtext = []
    for word in text.split():
        word = word.strip()
        while len(word) > maxchars:
            newtext.append(word[ : maxchars - 1] + '-')
            word = word[maxchars - 1 : ]
        newtext.append(word)
    return '\n'.join(newtext)


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
    plt.ylim(0, int(1.05 * ymax + 1))
    plt.xticks(xs, models, fontsize=32)
    plt.locator_params(axis='y', bins=4)
    plt.ylabel('sites out of %d (FDR %.2f)' % (nsites, fdr), fontsize=30)
    plt.legend(handles, labels, loc='upper right', numpoints=1, fontsize=34, borderaxespad=0)
    plt.savefig(plotfile, bbox_inches='tight')
    plt.clf()
    plt.close()



def SelectionViolinPlot(plotfile, ylabel, models, yvalues, symmetrizey, hlines=None, points=None, pointmarkercolor='or', usetex=True, legends=False, fixymin=None, fixymax=None, modelgroups=None):
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
    Alternatively, can be a single number -- in that case, we draw a horizontal
    line across the whole plot at that number.

    *points* : if not *None*, list of the same length as *models* with each entry
    a list giving the y-value for points to be placed for that *model*.

    *pointmarkercolor* : specifies marker and color of points in *points*. 
    Should either a length-two string giving marker and color for all points
    (such as *or* for circles, red) or a list of lists of the same length
    as *points* with each entry specifying the marker and color for that point.

    *usetex* : use LaTex formatting of strings?

    *legends* : Create legend(s) with names of points specified by *pointmarkercolor*?
    If it is not *False* or *None*, then *legends* should be a list. Each entry
    should be a 3-tuple *(markercolors, names, title)*. In this tuple,
    *markercolors* and *names* are lists of the same
    length, with *markercolors* being a list of marker / color (e.g. *or*
    for circles, red) and *names* being a list of the string corresponding
    to each marker / color; *title* is the title for the legend.

    *fixymin* : if not *None*, the y-minimum is fixed to this value.

    *fixymax* : if not *None*, the y-maximum is fixed to this value.

    *modelgroups* : do we "group" models on the x-axis? If so, set this to a list
    of the same length as *models* with each entry being the group to which
    that model is assigned. For instance, if *models* is 
    *['ExpCM', 'YNGKP', 'ExpCM', 'YNGKP']*,
    then *modelgroups* might be *['HA', 'HA', 'NP', 'NP']*. In this case,
    the two groups are indicated with a line and a label on the x-axis.
    Models in the same group must be consecutive. If any entry is *None*,
    the corresponding model is not assigned a group.
    """
    alpha = 0.55 # transparency for points
    markersize = 25 # size of points
    markerlw = 0.6 # line width for makers
    assert os.path.splitext(plotfile)[1].lower() == '.pdf', "plotfile %s does not end with extension '.pdf'"
    assert len(models) == len(yvalues) >= 1
    if modelgroups:
        assert len(modelgroups) == len(models), "modelgroups is not the same length as models"
        # make sure models in the same group are consecutive
        ngroups = 1
        previousgroup = modelgroups[0]
        for group in modelgroups[1 : ]:
            if group != previousgroup:
                ngroups += 1
                previousgroup = group
        assert ngroups == len(set(modelgroups)), "models in the same group must be consecutive in modelgroups. This is not the case:\n%s" % str(modelgroups)
    plt.rc('font', size=12)
    plt.rc('text', usetex=usetex)
    lmargin = 0.7
    tmargin = 0.1
    if modelgroups:
        bmargin = 0.6
    else:
        bmargin = 0.4
    if legends:
        perlegendwidth = 0.9
        rmargin = perlegendwidth * len(legends) + 0.03
    else:
        rmargin = 0.1
    (height, widthper) = (2.5, 1.5)
    violinwidth = 0.7
    totwidth = lmargin + rmargin
    if modelgroups:
        firstmodel = True
        withingroupspacing = violinwidth + 0.3 * (1.0 - violinwidth)
        xs = []
        for (imodel, igroup) in zip(models, modelgroups): 
            if not firstmodel and (igroup == lastgroup != None):
                xs.append(xs[-1] + withingroupspacing)
                totwidth += withingroupspacing * widthper
            elif firstmodel:
                firstmodel = False
                xs.append(0)
                totwidth += widthper
            else:
                xs.append(xs[-1] + 1)
                totwidth += widthper
            lastgroup = igroup
    else:
        xs = [x for x in range(len(models))]
        totwidth += widthper * len(models)
    totheight = height + tmargin + bmargin
    plt.figure(figsize=(totwidth, totheight))
    plt.axes([lmargin / totwidth, bmargin / totheight, 1.0 - (lmargin + rmargin) / totwidth, 1.0 - (tmargin + bmargin) / totheight])
    plt.ylabel(ylabel, fontsize=15)
    plt.violinplot(yvalues, xs, widths=violinwidth, showextrema=False)
    xmargin = 0.2 * violinwidth / 2.0
    xmin = xs[0] - violinwidth / 2.0 - xmargin
    xmax = xs[-1] + violinwidth / 2.0 + xmargin
    plt.xlim(xmin, xmax)
    if isinstance(hlines, (int, float)):
        plt.hlines(hlines, xmin, xmax, colors='b', linewidths=1, linestyles='dotted')
    elif hlines:
        assert len(hlines) == len(models)
        line_ys = []
        line_xmins = []
        line_xmaxs = []
        for (i, ix) in enumerate(xs):
            line_ys += hlines[i]
            line_xmins += [ix - violinwidth / 2.0] * len(hlines[i])
            line_xmaxs += [ix + violinwidth / 2.0] * len(hlines[i])
        plt.hlines(line_ys, line_xmins, line_xmaxs, colors='b', linewidths=1, linestyles='dotted')
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
        for (i, ix) in enumerate(xs):
            (model_xs, model_ys) = SmartJitter(points[i], yspace=(ymax - ymin) / 25., xspace=0.08, xcenter=ix)
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
    plt.xticks(xs, models, fontsize=15)
    if legends:
        legendx = 1.0 - rmargin / float(totwidth)
        legendfracwidth = perlegendwidth / float(totwidth)
        legendtop = 1.0 - tmargin / float(totheight)
        for (markercolors, legendnames, legendtitle) in legends:
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
            if legendtitle:
                legendtitle = ('\\bf{%s\nsites}' % SplitText(legendtitle, maxchars=6)).replace('\n', '}\n\\bf{')
            else:
                legendtitle = '\\bf{sites}'
            legend = plt.legend(sortedhandles, sortedlegendnames, bbox_to_anchor=(legendx, 0, legendfracwidth, legendtop), bbox_transform=plt.gcf().transFigure, fontsize=13, numpoints=1, title=legendtitle, markerscale=0.25, handlelength=0.7, handletextpad=0.25, borderaxespad=0, labelspacing=0.2)
            plt.gca().add_artist(legend)
            legendx += legendfracwidth
            plt.setp(legend.get_title(), fontsize=13)
    if modelgroups:
        for group in set(modelgroups):
            if not group:
                continue
            start_i = min([i for (i, g) in enumerate(modelgroups) if g == group])
            end_i = max([i for (i, g) in enumerate(modelgroups) if g == group])
            start_x = (xmargin + xs[start_i]) / (xmax - xmin) # axes coordinates
            end_x = (xmargin + xs[end_i] + violinwidth) / (xmax - xmin) # axes coordinates
            line_y = -0.12 # in axes coordinates
            cap_height = 0.03
            line = plt.Line2D([start_x, end_x], [line_y, line_y], transform=plt.gca().transAxes, color='black', linewidth=1.5, solid_capstyle='butt')
            line.set_clip_on(False)
            plt.gca().add_line(line)
            for x in [start_x, end_x]: # caps on end of lines
                line = plt.Line2D([x, x], [line_y + cap_height, line_y - cap_height], transform=plt.gca().transAxes, color='black', linewidth=1.5, solid_capstyle='butt')
                line.set_clip_on(False)
                plt.gca().add_line(line)
            plt.text((start_x + end_x) / 2.0, line_y - 0.04, group, transform=plt.gca().transAxes, horizontalalignment='center', verticalalignment='top', fontsize=15)
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
    if not ys:
        return ([], [])
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
