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
from scipy.signal import savgol_filter

parser = argparse.ArgumentParser()
parser.add_argument('ddir', default='', type=str, help='datadir directory')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='add verbosity')
parser.add_argument('-n', "--ncol", type=int, default=2, help="amount of columns in legend")
parser.add_argument('-l', "--light", action='store_true', default=False, help="Light color scheme")
parser.add_argument('-a', "--absolute", action='store_false', default=True, help="abs of helicity plot")
parser.add_argument('-e', '--error', action='store_true', default=False, help='helicity error')
args = parser.parse_args()
figdir = 'figures'

k0 = 80

if args.ddir[-1] == '/':
    args.ddir = args.ddir[:-1]

def round_exp(num, prec=0):
    formatstr = '%.{0}e'.format(prec)
    newnum = float( formatstr % (num))
    return newnum

def check_dirname_grid(dd):
    if dd.split('_')[0] in ['visc', 'prandtl', 'slope','nonhel','visc1024', 'prandtl1024']:
        return True
    else:
        return False

if args.verbose:
    for arg in vars(args):
        print('{:<9s}: {:<12s}'.format(arg, str(getattr(args, arg))))

def check_figdir(directory):
    if not path.exists(directory):
        mkdir(directory)
        if path.isdir(directory):
            return
    else:
        return
   
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

def setup_specax(ax):
    if args.absolute:
        ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(1/k0,dim.nxgrid//2/k0)
    ax.set_xlabel('$k/k_0$')# mode')
    if args.absolute: 
        ylims = (1e-15, 1e-7)
        ylabl = r'$\mid \mathcal{H}_k\mid$'
    else:
        ylims = (-2,2)
        ylabl = r'$\mathcal{H}_k$'
    if args.error:
        ylims = (1e-4, 1)
        ylabl = r'$|\mathcal{H}_k|k/E_k$'
    ax.set_ylabel(ylabl)
    #ylims = (round_exp(powerb[i+1,1]), round_exp( 1.5*powerb[i,:].max()))
    #ylims = (2e-10, 8e-5)
    if not args.absolute:
        ax.text(1.1/k0, 1.5, r'$\times 10^{-8}$')
    ax.set_ylim(ylims)
    if args.verbose:
        print('ylim : %.1e %.1e ' % (ax.get_ylim()))

def get_plottimes(tb):
    ptimes = []
    if args.ddir.startswith('delta'):
        pptimes = [0,.1, 1, 10,200]
    elif args.ddir.startswith('nonhel') and tb[-1] > 300:
        pptimes = [0,1,5,50,300]
    else:
        #pptimes = [1,10,50,300]
        pptimes = [0,1,10,100]
        #pptimes = [ round(x) for x in np.logspace(0,np.log10(tb[-1]), 4)]
    for t in pptimes:
        ti = search_indx(tb, t, eps= 0.05)
        if ti is not None:
            ptimes.append(ti)
    if args.verbose: 
        print('times plotted: ',pptimes)
        print('plot indices: ',ptimes)
    return ptimes

def get_string(ddir):
    import fractions as F
    #if all(dd.startswith('visc') for dd in dirs):
    if ddir.startswith('visc'):
        #strs = [r'$\nu_3 = %s$', r'$\nu_3 = %s$', r'$\nu = %s$', r'$\nu = %s$', r'$\nu = %s$', r'$\nu = %s$']
        if 'hyper' in ddir.split('_'):
            nu = '_3'
        else:
            nu = ''
        string = r'$\nu%s = %s$' % (nu, to_times(ddir.split('_')[-1]))
        #strings = [s % to_times(dd.split('_')[-1]) for s,dd in zip(strs,dirs)]
    elif ddir.startswith('prandtl'):
        string = 'Pr = %i'% int(float(ddir.split('_')[-1]))
    #    pass
    elif ddir.startswith('slope'):
        exp = ddir.split('_')[-1]
        if exp == '0.5':
            f1 = F.Fraction(exp)
            string = r'$E\sim k^{%i/%i}$' % (f1.numerator, f1.denominator)
        else:
            string = r'$E\sim k^{%s}$' % exp
    elif ddir.startswith('nonhel'):
        string = '$k_{max}=%s$' % ddir.split('_')[-1][1:]
    if args.verbose: 
        print(string)
    return string

def read_dim_krms(ddir):
    dim = pc.read_dim(datadir=ddir)
    krms = np.loadtxt(path.join(ddir, 'power_krms.dat')).flatten()[:dim.nxgrid//2]
    return dim, krms

def make_specplot(fig, ax, powerarr, tarr):
    plottimes = get_plottimes(tarr)
    clrindx = iter(np.linspace(0,1,len(plottimes)))
    kmax = []
    maxx = []
    for p, ph in enumerate(powerarr):
        if p in plottimes:
            if args.absolute:
                hel_plot = savgol_filter(np.abs(ph), 5, 1)[1:]
            else:
                hel_plot = savgol_filter(1e8*ph, 9,2)[1:]
            if args.error:
                hel_plot =(powerb[p,1:]/ (hel_plot*krms[1:]))**-1
            s = 't = %.0f'
            #s = 't = %.0f'
            curr_clr = cscheme(next(clrindx))
            print(hel_plot[0])
            ax.plot(krms[1:]/k0, hel_plot, label=s % tarr[p], color=curr_clr, linewidth=1.5)
            if args.verbose: print('t= %.1f, color: %.2f, %.2f, %.2f, %.2f' % (tarr[p], *curr_clr))

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
import matplotlib.patches as patches

cscheme, clr1, clr2, clr3 = select_colors()

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height

# Setup Figure and Axes
fwidth = 0.45
 
try:
    powerhel =  np.loadtxt(join(args.ddir, 'helspec_tot.dat'))
    th = np.loadtxt(join(args.ddir, 'thel_tot.dat'))
    if args.error:
        powerb = np.loadtxt(join(args.ddir,'powertot.dat'))
        tb = np.loadtxt(join(args.ddir, 'ttot.dat'))
except FileNotFoundError:
    th, powerhel = pc.read_power('powerhel_mag.dat', datadir=args.ddir)
    if args.error:
        tb, powerb = pc.read_power('power_mag.dat', datadir=args.ddir)
print('Shape of Power_hel: ', powerhel.shape)

dim, krms = read_dim_krms(args.ddir)
fig, ax  = newfig(fwidth)
if args.verbose:
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    print('spec figure width: %.1f, height: %.1f' % (width, height))

setup_specax(ax)
make_specplot(fig, ax, powerhel, th)

ax.legend(loc='lower left', ncol=args.ncol, frameon=False)
fig.tight_layout(pad=.3)

if args.ddir.endswith('/'):
    args.ddir = args.ddir[:-1]

plotname = '_helicity_mag_spec'
if args.absolute:
    plotname += '_abs'

if args.error:
    filename = args.ddir + '_hel_error'
else:
    filename = args.ddir + plotname
if args.verbose:
    print(filename)
savefig(fig, path.join(figdir, filename))
print('SUCCESS')
