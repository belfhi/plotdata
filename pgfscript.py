import numpy as np
import pencil as pc
from scipy.optimize import curve_fit
from helperfuncs import *
import matplotlib as mpl
import sys, os
mpl.use('pgf')

try:
    ddir = sys.argv[1]
except IndexError:
    print('must give a pencil code data directory:')
    print(' usage: python pgfplots.py "directory" ')

def figsize(scale):
    fig_width_pt = 523.5307                   # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "font.size": 10,
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt

# I make my own newfig and savefig functions
def newfig(width):
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename):
    plt.savefig('{}.pgf'.format(filename))
    plt.savefig('{}.pdf'.format(filename))

# read data
tb, powerb = pc.read_power('power_mag.dat', datadir=ddir)
print(tb.shape, tb[-1])

plottimes = []
for t in [0,.1, .2, .5, 1,2,5,10, 20, 50]:
    ti = search_indx(tb, t, eps= 0.01)
    plottimes.append(ti)
print(plottimes)

dim = pc.read_dim(datadir=ddir)
krms = np.loadtxt(os.path.join(ddir, 'power_krms.dat')).flatten()[:dim.nxgrid//2]
p0 = (-15.56, 2.20, 4.26)

# Simple plot
fig, ax  = newfig(0.6)
ax.set_ylim(1e-13, 1e-3)
ax.set_xlim(1,512)
ax.set_xlabel('k mode')
ax.set_ylabel(r'$E_{mag}$')
kmax = []
for p, pb in enumerate(powerb):
    xi = np.where( pb == pb.max())[0][0]; xi1 = xi - xi//3; xi2 = xi + xi//3
    po, pco = curve_fit(parabola, np.log10(krms[xi1:xi2]), np.log10(pb[xi1:xi2]), p0=p0)
    kmax.append(10**(po[1]))
    #print(xi, xi1, xi2)
    if p in plottimes:
        s = 't = %.1f' if tb[p] < 1.0 else 't = %.0f'
        ax.plot(pb, label=s % tb[p])
        ax.plot(krms[xi1], pb[xi1], "o", color='blue')
        ax.plot(krms[xi2], pb[xi2], "o", color='blue')
        ax.plot(krms[xi1:xi2], 10**parabola(np.log10(krms[xi1:xi2]), *po), color='blue')#, color=clrs[i])
ax.legend(loc = 'upper left')

savefig('nonhel_k_180_mag_spec')
