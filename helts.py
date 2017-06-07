#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import pencil as pc
from scipy.optimize import curve_fit
from helperfuncs import *
import matplotlib as mpl
mpl.use('pgf')
import sys
from mpl_toolkits.axes_grid1 import Grid
import argparse
from os import mkdir, path, listdir
from os.path import join, isdir, isfile
from fractions import Decimal, Fraction
from matplotlib.ticker import MultipleLocator

parser = argparse.ArgumentParser()
parser.add_argument('ddir', default='', type=str, help='datadir directory')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='add verbosity')
parser.add_argument('-l', '--light', action='store_true', default=False, help='use light color scheme')
parser.add_argument('-hel', '--helicity', action='store_true', default=False, help='do the plot of the helicity')
parser.add_argument('-a', '--a_rms', action='store_true', default=False, help='do the plot of <A>^2')

args = parser.parse_args()

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

def select_colors():
    if args.light:
        cscheme = plt.cm.Paired
        clr1 = (0.65, 0.81, 0.89, 1.00)
        clr2 = (0.93, 0.56, 0.28, 1.00)
        clr3 = (0.72, 0.61, 0.46, 1.00)
    else:
        cscheme = plt.cm.Dark2
        clr1 = (0.11, 0.62, 0.47, 1.00)
        clr2 = (0.61, 0.35, 0.65, 1.00)
        clr3 = (0.40, 0.40, 0.40, 1.00)
    return cscheme, clr1, clr2, clr3

def check_figdir(directory):
    if not path.exists(directory):
        mkdir(directory)
        if path.isdir(directory):
            return
    else:
        return
    
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
    "figure.figsize": figsize(0.95),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt

cscheme, clr1, clr2, clr3 = select_colors()

# check that there is a figure directory
figdir = 'figures'
check_figdir(figdir)

# Setup Figure and Axes
fwidth = 0.45
tser = pc.read_ts(datadir=args.ddir, quiet=not args.verbose)
par = pc.read_param(datadir=args.ddir, quiet=True)
kpeak = par.kpeak_aa
tau0 = (tser.brms[0]*kpeak)**-1
if args.verbose:
    print('tau = ', tau0)
#if args.peak:
#    args.tsplot = True
fig, ax  = newfig(fwidth)
ax.set_xscale('log')
ax.set_xlim(1e-2/tau0, 3e2/tau0)
#ax.set_yscale('log')
if args.helicity:
    ax.set_ylabel(r'$\mathcal{H}$')#/\mathcal{H}_0$')
    ax.plot(tser.t/tau0, np.abs(tser.ab_int)*10**5)#**2/tser.brms**2)#/tser.ab_int[0]))
    twx = ax.twinx()
    twx.set_ylabel(r'$\mathcal{H}/\mathcal{H}_0$')
    twx.plot(tser.t/tau0, tser.ab_int/tser.ab_int[0])#**2/tser.brms**2)#/tser.ab_int[0]))
elif args.a_rms:
    try:
        A2 = tser.arms**2
    except AttributeError:
        print('no arms in timeseries object')
        sys.exit(1)
    ax.plot(tser.t/tau0, A2*1e5)
    ax.text(0.02, 0.98, r'$\times 10^{-5}$', transform=ax.transAxes, va='top')
    ax.set_ylabel(r'$\langle A^2\rangle$')
    ax.set_ylim(1.0, 2.2)
    loc = MultipleLocator(base=0.2)
    ax.yaxis.set_major_locator(loc)


#ax.set_yscale
#ax.plot(aran)
#ax.set_ylim(-.88, -.82)
ax.text(0.02, 0.98, r'$\times 10^{-5}$', va='top', transform=ax.transAxes)
ax.set_xlabel(r'time $t/\tau_0$')
#fig.suptitle('Total Helicity ')
fig.tight_layout()
if args.ddir.endswith('/'):
    args.ddir = args.ddir[:-1]
if args.helicity:
    fname = join(figdir, args.ddir + '_helicity_ts')
elif args.a_rms:
    fname = join(figdir,  args.ddir + '_arms_ts')
else:
    print('must give either helicity or arms')
    sys.exit(1)
if args.verbose:
    print(fname)
savefig(fig, fname)
