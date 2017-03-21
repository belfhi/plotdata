#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import matplotlib as mpl
import sys
from scipy.integrate import romb, simps
import argparse
mpl.use('pgf')

parser = argparse.ArgumentParser()
parser.add_argument('ddir', default='', type=str, help='datadir directory')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='add verbosity')
parser.add_argument('-hel', '--helical', action='store_true', default=False, help='plot the helical line too')
args = parser.parse_args()

dstart = args.ddir.split('_')[0]

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

def nu_string(dd):
    ds = dd.split('_')
    if 'hyper' in ds:
        return r'$\nu_3=%s$' % to_times(ds[-1])
    elif 'nonhel' in ds:
        return r'$k_{max}=%s$' % ds[1][1:]
    elif dd == 'helical':
        return 'helical'
    elif 'prandtl' in ds:
        return r'Pr $={%d}$' % float(ds[-1])
    else:
        return r'$\nu=%s$' % to_times(ds[-1])

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


dirs = [s for s in listdir('.') if isdir(s) and s.startswith(dstart) and not isfile(join(s,'NOPLOT'))]
if len(dirs) > 1:
    try:
        dirs = sorted(dirs, key=lambda s: float(s.split('_')[-1]), reverse=True)
    except ValueError:
        dirs = sorted(dirs)
if args.helical:
    dirs.insert(0, 'helical')
if args.verbose:
    print(dirs)
clrindx = iter(np.linspace(0,1,len(dirs)))

# Simple plot
fig, ax  = newfig(0.45, ratio=0.75)


for dd in dirs:   
    if dstart.startswith('delta'):
        tstop = 1.
    else:
        tstop=0
    dim = pc.read_dim(datadir=dd)
    krms = np.loadtxt(join(dd,'power_krms.dat')).flatten()[:dim.nxgrid//2]
    if args.verbose:
        print('dir: ', dd)
    try:
        t = np.loadtxt(join(dd,'ttot.dat'))
        powerb = np.loadtxt(join(dd,'powertot.dat'))
        if args.verbose:
            print('reading concatenated data files')
    except FileNotFoundError:
        t, powerb = pc.read_power('power_mag.dat', datadir=dd)
    tsind = search_indx(t, tstop)
    if args.verbose: 
        print('tstop: %.1f, index = %i' % (tstop, tsind))
    #if args.verbose:
    #    print('%.1f %.1f %.1f' %(t[0], t[1], t[-1]))
    kmax = 7 # arbitrary value
    dim = pc.read_dim(datadir=dd)
    emax = []
    for p,pb in enumerate(powerb):
        #if round(t[p],1) < tstop:
        #    continue
        a = simps(pb[1:kmax], krms[1:kmax])
        emax.append(a)
    emax1 = emax/emax[0]
    #tvals = len(t)-len(emax1)
    t = t[len(t)-len(emax1):]
    #ll = '$k^{%g}$' % float(dd[1:]))
    #ll = r'$\nu=%s$' % to_times(dd.split('_')[-1]))
    if args.verbose:
        print(t.shape, emax1.shape)
    clr = plt.cm.Dark2(next(clrindx))
    if args.verbose:
        print('clr: %.1f, %.1f, %.1f, %.1f' % clr)
    #ax.plot(t[tvals:], emax1, label=nu_string(dd), color=clr, linewidth=1.5)
    ax.plot(t[1:], emax1[1:], label=nu_string(dd), color=clr, linewidth=1.5)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(tstop,t[-1])
if not args.helical:
    ax.set_ylim(1e-1, 1e2)
#ax.set_yticks([10,20,50,100])
#ax.set_yticklabels(['$10$','$20$','$50$','$100$'])
if not args.helical:
    ax.legend(loc='lower left', ncol=2, frameon=False)
else:
    ax.legend(loc='upper left', ncol=1, fancybox=True, framealpha=0.5)
#fig.suptitle('Wavenumber of Maximum ')

ax.set_xlabel('time')
ax.set_ylabel(r'$E_{k\leq k_7}/E_0$')
fig.tight_layout(pad=0.3)
if args.helical:
    filename = dstart + '_energy_increase_hel'
else:
    filename = dstart + '_energy_increase'
if args.verbose:
    print(filename)
savefig(join('figures',filename))
print('SUCCESS')
