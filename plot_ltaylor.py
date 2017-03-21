#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import pencil as pc
from scipy.optimize import curve_fit
from helperfuncs import *
import matplotlib as mpl
mpl.use('pgf')
import sys
from os import listdir, path
from mpl_toolkits.axes_grid1 import Grid
import argparse
from os import mkdir, path, listdir
from os.path import join, isdir, isfile

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='add verbosity')
parser.add_argument('ddir', default='', type=str, help='datadir directory')
parser.add_argument('choice', default='', type=str, choices=('taylor', 'energy'), help='choice of plotting')

args = parser.parse_args()
figdir = 'figures'

dstart = args.ddir.split('_')[0]
dirs = [l for l in listdir('.') if l.startswith(dstart) and isdir(l)]
try:
    dirs = sorted(dirs, key=lambda a: float(a.split('_')[-1]))
except ValueError:
    dirs = sorted(dirs, key=lambda a: float(a.split('_')[-1][1:]))
dirs.append('helical')
if args.verbose: print(dirs)

#if args.choice == 'energy':
#    print('energy plot not yet implemented, choose taylor instead')
#    exit(1)

def figsize(scale, ratio=None):
    #fig_width_pt = 418.25555
    fig_width_pt = 523.5307              # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    if not ratio:
        golden_mean = (np.sqrt(5.0)-1.0)/2.0   #golden_mean          
        ratio = golden_mean
    fig_height = fig_width*ratio             # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

# I make my own newfig and savefig functions
def newfig(width):
#    plt.clf()
    fig = plt.figure(figsize=figsize(width, ratio=0.75))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(fig, filename):
    fig.savefig('{}.pgf'.format(filename))
    fig.savefig('{}.pdf'.format(filename))

def make_tsplot(fig, ax):
    if args.choice == 'taylor':
        try:
            plotarr = 2*np.pi*tser.brms/tser.jrms
        except AttributeError:
            return
    elif args.choice == 'energy':
        plotarr = tser.brms**2
    ax.plot(tser.t, plotarr, label=dd.split('_')[-1], linewidth=1.5, 
            color=plt.cm.Paired(next(clrindx)))
    ti = search_indx(tser.t, 1., eps=.05)
    ti2 = search_indx(tser.t, 0.8*tser.t[-1], eps=.05)
    if args.verbose: print('index from which fitting begins: ',ti)
    if ti is not None:
        po1, pco1 = curve_fit(lambda x, a, b: powerlaw(x, a, b, x0=0), tser.t[ti:], plotarr[ti:])
        #ax.plot(tser.t[ti:ti2], powerlaw(tser.t[ti:ti2],1.2*po1[0],po1[1]), linewidth=1.5, 
        #        label=r'fit: $B\sim t^{%.1f}$' % po1[1], color=(0.72, 0.61, 0.46, 1.00))
        if args.verbose:
            print(dd, 'L ~ t^{1:.2f}'.format(*po1))

def setup_tsplot(ax):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('time')
    ax.set_ylabel(r'$L_{\textrm{T}}$')

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all tex
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "font.size": 10,
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.figsize": figsize(0.95),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt


# Setup Figure and Axes

fig, ax  = newfig(0.45)
setup_tsplot(ax)#, xlim=(.1,tser.t[-1]))
tmax = 0.
clrindx = iter(np.linspace(0,1,len(dirs)))

for dd in dirs:
    if args.verbose: print(dd)
    tser = pc.read_ts(datadir=dd, quiet=not args.verbose)
    make_tsplot(fig, ax)
    if tser.t[-1] > tmax:
        ax.set_xlim(.1, tser.t[-1])
    ax.legend(loc='upper left', frameon=False)

filename = dstart + '_taylor'
fig.tight_layout(pad=0.3)
if args.verbose:
    print(filename)
savefig(fig, path.join(figdir, filename))
print('SUCCESS')
