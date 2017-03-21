#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import pencil as pc
from scipy.optimize import curve_fit
from scipy.integrate import romb, simps
from helperfuncs import *
import matplotlib as mpl
mpl.use('pgf')
import sys
import argparse
from os import mkdir, path, listdir
from os.path import join, isdir, isfile

parser = argparse.ArgumentParser()
parser.add_argument('ddir', choices=['visc', 'nonhel', 'prandtl', 'delta', 'slope'],
                    type=str, help='directory from which data is read')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='add verbosity')
parser.add_argument('-t', '--tplot', action='store', default=10, type=int, help='specify plotting time')
parser.add_argument('-e', '--energy', action='store_true', default=False, help='do energy plot, not length scale plot')
parser.add_argument('-hel', '--helical', action='store_true', default=False, help='add a data point for helicity run')
parser.add_argument('-twin', action='store_true', default=False, help='add twin axes')

args = parser.parse_args()
figdir = 'figures'

fhyper = True
fnu = True

if args.twin:
    fnameadd = 'twin'
elif args.energy:
    fnameadd = 'energy'
else:
    fnameadd = 'lint'

def sortdirs(dd, flt=True):
    if not args.ddir == 'nonhel':
        s = dd.split('_')[-1]
    else:
        s = dd.split('_')[-1][1:]
    if flt:
        return float(s)
    else:
        return s

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

def to_times(s):
    s1, s2 = s.split('e')
    return r'{0}{{\times}} 10^{{{1}}}'.format(s1, s2)

def get_labels(hyper3=False, energy=False):
    if args.twin:
        if energy:
            ll = r'$E_{k\leq 7}/E_0$'
            clr = (0.65, 0.81, 0.89, 1.00)
        else:
            ll = r'$L_{\textrm{corr}}$'
            clr = (0.93, 0.56, 0.28, 1.00)
    elif args.ddir == 'visc':
        if hyper3:
            ll = r'$\nu_{\textrm{hyper3}}$'
            clr = (0.65, 0.81, 0.89, 1.00)
        else:
            ll = r'$nu_{\textrm{const}}$'
            clr = (0.93, 0.56, 0.28, 1.00)
    else:
        ll = ''
        clr = (0.93, 0.56, 0.28, 1.00)
    return ll, clr

dirs = sorted([l for l in listdir('.') if l.startswith(args.ddir) and isdir(l) and not isfile(join(l,'NOPLOT'))], key=sortdirs)
if args.helical:
    dirs = ['helical'] + dirs
if args.verbose:
    print(dirs)

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

fwidth = 0.45
fig, ax = newfig(fwidth)
if args.twin:
    twinax = ax.twinx()
    twinax.set_yscale('log')
    twinax.set_ylabel(r'$E_{k\leq 7}/E_0$')

first_dir = True

# lists to store the data
Ek7  = []
kmax = []
#special lists for visc different colors:
ehyper = []

for i,dd in enumerate(dirs):
    if args.verbose: 
        print('current dir: ', dd)

    # starting value for fitting process
    p0 = (-12.92, 1.94, 3.46)  #(-15.56, 2.20, 4.26)
    
    # read in parameter and krms file
    dim = pc.read_dim(datadir=dd)
    krms = np.loadtxt(join(dd, 'power_krms.dat')).flatten()[:dim.nxgrid//2]
    xi = dim.nxgrid//2
    par2 = pc.read_param(quiet=True, datadir=dd, param2=True)

    # reading the spectrum data
    try:
        tb, powerb = (np.loadtxt(join(dd, 'ttot.dat')), np.loadtxt(join(dd, 'powertot.dat')))
    except FileNotFoundError:
        tb, powerb = pc.read_power('power_mag.dat', datadir=dd) 
    if args.verbose: 
        print('n_spec = %g; k_Ny = %g' % powerb.shape)

    # look for plotting index
    indx = search_indx(tb, args.tplot)
    if args.verbose: 
        print('selected plotting index: ',indx)

    # dont do analysis for first helical spectra 
    # should be removed, when proper restart data is there
    if dd == 'helical':
        tstop = 0.5
    else:
        tstop = 0
    if args.verbose: 
        print('t_stop = ', tstop)

    # start loop over data
    first = True
    for p,pb in enumerate(powerb):
        if args.energy or args.twin:
            km = 7
            if round(tb[p],1) >= tstop:
                if first:
                    E0 = simps(pb[1:km], krms[1:km])
                    first = False
                    print('first energy: ', E0)
 
                elif p == indx:
                    etemp = simps(pb[1:km], krms[1:km])/E0
                    if args.ddir == 'visc' and any( x in dd.split('_') for x in ['helical', 'hyper']):
                        ehyper.append(etemp)
                    else:
                        Ek7.append(etemp)

        if not args.energy or args.twin:
            if round(tb[p],1) >= tstop:
                xi = np.where( pb == pb[:xi+1].max())[0][0]; xi1 = xi - xi//3; xi2 = xi + xi//2
                try:
                    po, pco = curve_fit(parabola, np.log10(krms[xi1:xi2]), np.log10(pb[xi1:xi2]), p0=p0)
                    if p == 0 and args.verbose and first_dir:
                        print('fitting parameters:  \n',po)
                except RuntimeError:
                    print('fitting not successfull')
                    print(xi, xi1, xi2)
                    sys.exit(1)
                p0 = po
                if p == 0:
                    km0 = 10**(-po[1])
                if p == indx:
                    kmax.append(10**(-po[1]))
            else:
                continue
    first_dir = False

mkr = 6
xvals =[i+1 for i in range(len(dirs))]
if args.verbose:
    print('xvals: ', xvals, ' len :', len(xvals))
if args.energy or args.twin:
    if args.twin:
        ll, clr = get_labels(energy=True)
        l1 = twinax.plot(xvals, Ek7, 'o', color=clr, markersize=mkr, label=ll)
    else:
        if any('hyper' in dd.split('_') for dd in dirs):
            ll, clr = get_labels(hyper3=True)
            if args.helical:
                l1 = ax.plot([1,2,3], ehyper, 'o', color=clr, markersize=mkr, label=ll)
                ll, clr = get_labels(hyper3=False)
                l2 = ax.plot([4,5,6,7], Ek7, 'o', color=clr, markersize=mkr, label=ll)
            else:
                l1 = ax.plot([1,2], ehyper)   
                ll, clr = get_labels(hyper3=False)
                l2 = ax.plot([3,4,5,6], Ek7, 'o', color=clr, markersize=mkr, label=ll)
        else:
            l1 = ax.plot(xvals, Ek7, 'o', color=clr, markersize=mkr, label=ll)
    print('plotting energy')

if not args.energy or args.twin:
    ll, clr = get_labels(energy=False)
    l2 = ax.plot(xvals, kmax, 'o', color=clr, markersize=mkr, label=ll)
    print('plotting length scale')

if args.verbose:
    print('len(E_k7) = %g; len(tb) = %g '% (len(Ek7), *tb.shape))

ylim = ax.get_ylim()
#ax.set_xlim(0, len(dirs)+1)
ax.set_xlim(0.75, len(dirs)+0.25)
if args.helical:
    ax.set_yscale('log')
    #print('ylims: {:.1f} {:.1f}'.format( *ax.get_ylim()))
    #ax.set_ylim(0.5*ylim[0], 1.1*ylim[1])
else:
    ax.set_ylim(0.9*ylim[0], 1.1*ylim[1])
if args.verbose:
    print('ylims: {:.1f} {:.1f}'.format( *ax.get_ylim()))
ax.set_xticks(range(1,len(dirs)+1))

if args.ddir == 'visc' and not args.energy:
    txtloc = (0.8,0.05)
elif args.ddir == 'visc' and args.energy:
    txtloc = (0.05, 0.05)
elif args.twin:
    txtloc = (0.05, 0.9)
else:
    txtloc = (0.8,0.9)

ax.text(*txtloc, 't=%.0f'% args.tplot, transform=ax.transAxes, horizontalalignment='left')

if args.twin:
    ax.set_ylabel(r'$L_{\textrm{corr}}$')
else:
    if args.energy:
        ax.set_ylabel(r'$E_{{k\leq k_{{{0:d}}}}}/E_0$'.format(km))
    else:
        ax.set_ylabel(r'$L_{\textrm{corr}}$')

if args.ddir == 'prandtl':
    #ax.set_xticklabels(['$%s$' % to_times(sortdirs(s, flt=False)) for s in dirs])
    if args.helical:
        xtlabls = ['%.0f' % float(s.split('_')[-1])  for s in dirs[1:]]
        xtlabls.insert(0,r'$\mathcal{H}$')
    else:
        xtlabls = ['%.0f' % float(s.split('_')[-1])  for s in dirs]
    ax.set_xlabel('Prandtl number')
elif args.ddir == 'visc':
    #ax.set_xticklabels(['$%s$' % to_times(sortdirs(s, flt=False)) for s in dirs], rotation=45)
    if args.helical:
        xtlabls = ['$%s$' % to_times(sortdirs(s, flt=False)) for s in dirs[1:]]
        xtlabls.insert(0, r'$\mathcal{H}$')
    else:
        xtlabls = ['$%s$' % to_times(sortdirs(s, flt=False)) for s in dirs]
    ax.set_xlabel(r'viscosity parameter $\nu$')
else:
    xtlabls =['%.0e' % sortdirs(s) for s in dirs]

ax.set_xticklabels(xtlabls)
if args.twin or args.ddir == 'visc':
    lns = l1 + l2
    lbls = [ln.get_label() for ln in lns]
    ax.legend(lns, lbls, loc='upper right', numpoints=1)
#else:
#    ax.legend(loc='lower left', numpoints=1)

if args.helical:
    fname = '%s_t%.0f_%s_comparison_hel' % (args.ddir, args.tplot,fnameadd )
else:
    fname = '%s_t%.0f_%s_comparison' % (args.ddir, args.tplot,fnameadd )
fname = join(figdir, fname)
if args.verbose:
    print(fname)
#ax.autoscale(tight=True) 
fig.tight_layout(pad=0.3)
savefig(fig, fname)
print('SUCCESS')
