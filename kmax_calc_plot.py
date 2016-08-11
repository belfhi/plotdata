#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import matplotlib as mpl
import sys
import argparse
mpl.use('pgf')

parser = argparse.ArgumentParser()
parser.add_argument('ddir', default='', type=str, help='datadir directory')
#parser.add_argument('-s', '--spectra',  action='store_true', default=False, help="plot magnetic spectra")
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='add verbosity')
args = parser.parse_args()

def figsize(scale, ratio=None):
    fig_width_pt = 523.5307                         # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    if  not ratio:    
        golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
        fig_height = fig_width*golden_mean 
    else:
        fig_height = fig_width*ratio
        # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "font.size": 10,
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

def newfig(width, ratio=None):
    plt.clf()
    fig = plt.figure(figsize=figsize(width, ratio))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename):
    plt.savefig('{}.pgf'.format(filename))
    plt.savefig('{}.pdf'.format(filename))

def to_times(s):
    s1, s2 = s.split('e')
    return r'{0}\times 10^{{{1}}}'.format(s1, s2)

def linear(x, a, b):
    return a*x + b

def powerlaw(x, a, b):
    return a*x**b

def parabola(x, a, b, c):
    return a*(x-b)**2 + c

import pencil as pc
import numpy as np
from helperfuncs import search_indx
from os import listdir, path
from os.path import isdir, isfile, join
from scipy.optimize import curve_fit


dirs = [d for d in listdir('.') if d.startswith(args.ddir[:4]) and isdir(d) and not isfile(join(d,'NOPLOT'))]
dirs = sorted(dirs, key=lambda a: float(a.split('_')[-1]), reverse=True)
#dirs =sorted( [s for s in listdir('.') if (isdir(s) and s.startswith(dstart) and not isfile(join(s,'.noscale')))],
             #key=lambda s: float(s.split('_')[-1]))
dstart=dirs[0].split('_')[0]
if args.verbose:
    print(dirs)
clrindx = iter(np.linspace(0,1,len(dirs)))

dim = pc.read_dim(datadir=args.ddir)
krms = np.loadtxt(join(args.ddir, 'power_krms.dat')).flatten()[:dim.nxgrid//2]
if args.verbose:
    print('krms shape: ',krms.shape)


# Simple plot
fig, ax  = newfig(0.45, ratio=0.75)
ax.set_xscale('log')
first_dir = True
for dd in dirs:   
    p0 = (-12.92758524, 1.94666781, 3.4643292)  #(-15.56, 2.20, 4.26)
    if args.verbose:
        print('dir: ', dd)
    kmax = []
    dim = pc.read_dim(datadir=dd)
    xi = 256
    try:
        t = np.loadtxt(join(dd, 'tb_res.dat'))
        powerb = np.loadtxt(join(dd, 'powerb_res.dat'))
    except FileNotFoundError:
        t, powerb = pc.read_power('power_mag.dat', datadir=dd)
    if args.verbose: 
        print(t.shape, powerb.shape)
    for p,pb in enumerate(powerb):
        xi = np.where( pb == pb[:xi+1].max())[0][0]; xi1 = xi - xi//3; xi2 = xi + xi//2
        try:
            po, pco = curve_fit(parabola, np.log10(krms[xi1:xi2]), np.log10(pb[xi1:xi2]), p0=p0)
            if p == 0 and args.verbose and first_dir:
                print(po)
        except RuntimeError:
            print(xi, xi1, xi2)
        p0 = po
        kmax.append(1./(10**po[1]))
    if args.verbose:
        print('len(kmax, t): ',len(kmax), t.shape)
    if dd.startswith('visc'):
        if dd.split('_')[2] == 'hyper':
            hyp3 = '_3'
        else:
            hyp3 = ' '
        s = r'$\nu%s = %s$' % (hyp3, to_times(dd.split('_')[-1]))
    else:
        s = r'Pr $=%i$' % float(dd.split('_')[-1])
    ax.plot(t[1:], kmax[1:], label=s  , linewidth=1.5, color=plt.cm.Paired(next(clrindx)))
    first_dir = False
ax.set_yscale('log')
ax.set_xlim(.1,100)
if dstart.startswith('visc_'):
    ax.set_ylim(1e-2, 2e-1)
else:
    ax.set_ylim(1e-2, 1e-1)
#ax.set_yticks([10,20,50,100])
#ax.set_yticklabels(['$10$','$20$','$50$','$100$'])
ax.legend(loc='upper left', ncol=1, frameon=False)
#fig.suptitle('Wavenumber of Maximum ')

ax.set_xlabel('time')
ax.set_ylabel(r'$L_{\textrm{int}}$')
fig.tight_layout(pad=0.3)
filename = dstart +'_scale_evolution'
if args.verbose: print(filename)
savefig(filename)
print('SUCCESS')
