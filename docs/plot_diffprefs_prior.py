"""Plots differential preferences prior."""


import math
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def F(x, c1, c2):
    return -c2 * math.log(1.0 + c1 * x**2)

plotfile = 'diffprefs_prior_plot.png'

plt.rc('text', usetex=True)
plt.rc('font', size=13)
xmax = 1
x = 0
xs = [x]
while x <= xmax:
    xs.append(x)
    x += 0.001
labels = []
handles = []
plt.figure(figsize=(9, 5))
lmargin = 0.11
rmargin = 0.32
bmargin = 0.13
tmargin = 0.07
plt.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
for (c1, color) in zip([50, 150, 500], ['b', 'r', 'g', 'y']):
    for (c2, ls) in zip([0.25, 0.5, 1], ['-', ':', '--']):
        ys = [F(x, c1, c2) for x in xs]
        label = '$C_1 = %d$, $C_2 = %.2f$' % (c1, c2)
        labels.append(label)
        handle = plt.plot(xs, ys, label=label, ls=ls, color=color, lw=2)
        handles.append(handle[0])
plt.ylabel('log likelihood: $-C_2 \\times \log \left(1 + C_1 \\times \left(\Delta\pi_{r,a}\\right)^2\\right)$', size=17)
plt.xlabel('$|\Delta\pi_{r,a}| = |\hat{\pi}_{r,a} - \left(\pi_{r,a}\\right)^{\\beta}|$', size=17)
plt.title('Prior over differential preferences', size=17)
plt.figlegend(handles, labels, loc='center right')
plt.show()

plt.savefig(plotfile)
