"""Plots differential preferences prior."""


import math
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

def F(x, c):
    return math.log(1.0 / (1.0 + c * x))

plotfile = 'diffprefs_prior_plot.png'

plt.rc('text', usetex=True)
plt.rc('font', size=13)
xmax = 2
x = 0
xs = [x]
while x <= xmax:
    xs.append(x)
    x += 0.001
labels = []
handles = []
plt.figure(figsize=(8, 5))
lmargin = 0.12
rmargin = 0.25
bmargin = 0.13
tmargin = 0.07
plt.axes([lmargin, bmargin, 1.0 - lmargin - rmargin, 1.0 - tmargin - bmargin])
for (c, color, ls) in zip([1, 10, 100, 1000], ['b', 'r', 'g', 'c'], ['-', ':', '--', '-.']):
        ys = [F(x, c) for x in xs]
        label = '$C = %d$' % c
        labels.append(label)
        handle = plt.plot(xs, ys, label=label, ls=ls, color=color, lw=2)
        handles.append(handle[0])
plt.ylabel('log likelihood: $\log \left(\\frac{1}{1 + C \sum_a \left(\Delta\pi_{r,a}\\right)^2}\\right)$', size=17)
plt.xlabel('$\sum_a \left(\Delta\pi_{r,a}\\right)^2$', size=17)
plt.title('Inverse-square prior over differential preferences', size=17)
plt.figlegend(handles, labels, loc='center right')
plt.show()

plt.savefig(plotfile)
