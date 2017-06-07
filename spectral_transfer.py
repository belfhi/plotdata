i#!/usr/bin/env python3
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

# check that there is a figure directory
figdir = 'figures'
check_figdir(figdir)

# read data from T_kpq.np
T_kpq = np.load(join(args.ddir,'T_kpq.npy'))
print(T_kpq.shape)

i = 8
j = 1
ex = 150
sigma = 2
k0=80

xlbl = ['$p/p_0$','$q/q_0$', '$k/k_0$']
ylbl = ['$q/q_0$','$k/k_0$', '$p/p_0$']
ax2.text(0.9,0.05, '{0} = {1!s}'.format(xlbl[j-1], i/k0), transform=ax2.transAxes, ha='right')
Tk = T_kpq[i,:,:]
Tp = T_kpq[:,i,:]
Tq = T_kpq[:,:,i]

Tk_smooth =0.5*( Tk[::2,::2] + Tk[1::2,1::2])
Tp_smooth =0.5*( Tp[::2,::2] + Tp[1::2,1::2])
Tq_smooth =0.5*( Tq[::2,::2] + Tq[1::2,1::2])

T_list = [Tk, Tp, Tq]
T_smooth = [Tk_smooth, Tp_smooth, Tq_smooth]
gauss_image = gaussian_filter(T_list[j], sigma)
#ax2.imshow(np.real(T_kpq[256,:,:]), cmap=plt.cm.hot)
im = ax2.imshow(gauss_image[:ex,:ex], extent=[0,ex/k0,ex/k0,0])
#divider = make_axes_locatable(ax2)
#ax2.set_xlim(0,100/k0)
#ax2.set_ylim(0,100/k0)
#cax = divider.append_axes("right", size="5%", pad=0.05)
ax2.set_xlabel(xlbl[j])
ax2.set_ylabel(ylbl[j])
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')
# Setup Figure and Axes
fwidth = 0.45
#tser = pc.read_ts(datadir=args.ddir, quiet=not args.verbose)
par = pc.read_param(datadir=args.ddir, quiet=True)
kpeak = par.kpeak_aa
tau0 = (tser.brms[0]*kpeak)**-1
if args.verbose:
    print('tau = ', tau0)
#if args.peak:
#    args.tsplot = True
fig, ax  = newfig(fwidth)

ax.text(0.92, 0.08, r'$\times 10^{-5}$', va='top', transform=ax.transAxes)
ax.set_xlabel(r'time $t/\tau_0$')
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
