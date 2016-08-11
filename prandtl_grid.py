#!/usr/bin/env python3
import numpy as np
import pencil as pc
from scipy.optimize import curve_fit
from helperfuncs import *
import matplotlib as mpl
mpl.use('pgf')
import sys
from os import listdir, path
from mpl_toolkits.axes_grid1 import Grid

def fsize(scale):
    fig_width_pt = 523.5307              # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*2/3.  #golden_mean         # height in inches
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
    "figure.figsize": fsize(0.95),     # default fig size of 0.9 textwidth
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

def to_times(s):
    s1, s2 = s.split('e')
    return r'{0}\times 10^{{{1}}}'.format(s1, s2)

def plot_grid_spectra(ax, i):
    #print(i)
    plottimes = []
    for t in [0,1,2,5,10, 20, 50, 100,200]:
        ti = search_indx(tb_ges[i], t, eps= 0.02)
        if ti is not None:
            plottimes.append(ti)
    #print('%s has %i lines' % (dirs[i],len(plottimes)))
    #print(plottimes)
    for p, pb in enumerate(pb_ges[i]):
        if p in plottimes:
            #print('plotting',dd, p)
            ax.loglog(krms, pb)
    ax.set_xlim(1,dim.nxgrid//2)
    ax.set_ylim(pb_ges[i][0,1], 1.5*pb_ges[i][0,:].max())

dirs = [d for d in listdir('.') if d.startswith('prandtl') and path.isdir(d)]
dirs = sorted(dirs, key=lambda a: float(a.split('_')[-1]))
print(dirs)

pb_ges = []
tb_ges = []
pr_numbers = {}
for dd in dirs:
    tb, powerb = pc.read_power('power_mag.dat', datadir=dd)
    pr_numeric = float(dd.split('_')[-1])
    #print(nu_numeric)
    pr_numbers[dd] = pr_numeric
    pb_ges.append(powerb)
    tb_ges.append(tb)

#print('tb_ges has : %i indices' % len(tb_ges))

dim = pc.read_dim(datadir=dirs[0])
krms = np.loadtxt(path.join(dirs[0],'power_krms.dat')).flatten()[:dim.nxgrid//2]

fig = plt.figure(figsize=fsize(.95))#figsize(width)

grid = Grid(fig, rect=[.08,.1,.9,.85], nrows_ncols=(2,2),
            axes_pad=0., label_mode='L',
            )
strings = [r'$\textrm{Pr} = %i$', r'$\textrm{Pr} = %i$', r'$\textrm{Pr} = %i$', r'$\textrm{Pr} = %i$']

for i,ax in enumerate(grid):
    plot_grid_spectra(ax, i)
    ax.title.set_visible(False)
    ax.text(1.5, 1e-4, strings[i] % int(float(dirs[i].split('_')[-1])))


fig.text(0.52, 0.01, 'k mode', ha='center')
fig.text(0.01, 0.5, r'$E_{mag}$', va='center', rotation='vertical')

savefig('grid_prandtl_spec')

