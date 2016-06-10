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
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('ddir', default='', type=str, help='datadir directory')
parser.add_argument('-s', '--spectra',  action='store_true', default=False, help="don't plot spectra")
parser.add_argument('-t', '--tsplot', action='store_true', default=False, help="plot timeseries too")
parser.add_argument('-g', '--gridplot', action='store_true', default=False, help='do the gridplot for multiple dirs')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='add verbosity')
parser.add_argument('-f', '--fit', action='store_false', default=True, help="don't fit a parabola")
parser.add_argument('-p', '--peak', action='store_true', default=False, help="make a plot of peak evolution")

args = parser.parse_args()

if args.verbose:
    for arg in vars(args):
        print(arg, getattr(args, arg))

if args.gridplot and not ( args.ddir.startswith('visc_') or args.ddir.startswith('prandtl')):
    print('gridplot only works for nu comparison and prandtl numbers')
    print('usage:')
    print('   python pgfscript.py -g visc_')
    print(' or:')
    print('   python pgfscript.py -g prandtl')
    sys.exit(1)
elif args.gridplot and (args.ddir.startswith('visc_') or args.ddir.startswith('prandtl')):
    dirs = [d for d in listdir('.') if (d.startswith(args.ddir[:4]) and path.isdir(d))]
    dirs = sorted(dirs, key=lambda a: float(a.split('_')[-1]))
    if args.verbose: print(dirs)

    
def read_grid_data(dirs):
    pb_ges = []
    tb_ges = []
    dd_visc = {}
    for dd in dirs:
        tb, powerb = pc.read_power('power_mag.dat', datadir=dd)
        nu_numeric = float(dd.split('_')[-1])
        #print(nu_numeric)
        dd_visc[dd] = nu_numeric
        pb_ges.append(powerb)
        tb_ges.append(tb)
    return pb_ges, tb_ges
    

def figsize(scale, ratio=None):
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
    return r'{0}\times 10^{{{1}}}'.format(s1, s2)

def setup_specax(ax):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(powerb[1,1], 1.5*powerb[0,:].max())
    ax.set_xlim(1,dim.nxgrid//2)
    ax.set_xlabel('k mode')
    ax.set_ylabel(r'$E_{mag}$')

def get_plottimes(tb):
    ptimes = []
    #for t in [0,.1,1,10,200]:
    for t in [0, 1, 10, 100]:
        ti = search_indx(tb, t, eps= 0.02)
        if ti is not None:
            ptimes.append(ti)
    if args.verbose: print(ptimes)
    return ptimes

def get_strings(dirs):
    if all(dd.startswith('visc') for dd in dirs):
        strs = [r'$\nu_3 = %s$', r'$\nu_3 = %s$', r'$\nu = %s$', r'$\nu = %s$', r'$\nu = %s$', r'$\nu = %s$']
        strings = [s % to_times(dd.split('_')[-1]) for s,dd in zip(strs,dirs)]
    elif all(dd.startswith('prandtl') for dd in dirs):
        strs = ['Pr = %i']*len(dirs)
        strings = [s % int(float((dd.split('_')[-1]))) for s,dd in zip(strs,dirs)]
    #    pass
    if args.verbose: 
        print(strings)
    return strings

def read_dim_krms(ddir):
    dim = pc.read_dim(datadir=ddir)
    krms = np.loadtxt(path.join(ddir, 'power_krms.dat')).flatten()[:dim.nxgrid//2]
    return dim, krms

def plot_grid_spectra(ax, i):
    plottimes = get_plottimes(tb_ges[i])
    clrindx = iter(np.linspace(0,1,len(plottimes)))
    for p, pb in enumerate(pb_ges[i]):
        if p in plottimes:
            ax.loglog(krms, pb, linewidth=1.5, color=plt.cm.Paired(next(clrindx)))
    ax.set_xlim(1,dim.nxgrid//2)
    ax.set_ylim(pb_ges[i][0,1], 1.5*pb_ges[i][0,:].max())

def make_specplot(fig, ax):
    plottimes = get_plottimes(tb)
    clrindx = iter(np.linspace(0,1,len(plottimes)))
    kmax = []
    maxx = []
    for p, pb in enumerate(powerb):
        if p == 0 and args.ddir.startswith('delta'):
            continue
        xi = np.where( pb == pb.max())[0][0]; xi1 = xi - xi//3; xi2 = xi + xi//3
        if args.fit:
            po, pco = curve_fit(parabola, np.log10(krms[xi1:xi2]), np.log10(pb[xi1:xi2]), p0=p0)
            kmax.append(10**(po[1]))
            maxx.append(max(pb))
        #print(xi, xi1, xi2)
        if p in plottimes:
            #s = 't = %.1f' if tb[p] < 1.0 else 't = %.0f'
            s = 't = %.0f'
            ax.plot(krms[1:], pb[1:], label=s % tb[p], color=plt.cm.Paired(next(clrindx)), linewidth=1.5)
            #ax.plot(krms[xi1], pb[xi1], "o", color='blue')
            #ax.plot(krms[xi2], pb[xi2], "o", color='blue')
            #ax.plot(krms[xi1:xi2], 10**parabola(np.log10(krms[xi1:xi2]), *po), color='blue')#, color=clrs[i])
    if args.peak:
        fig2, ax2 = newfig(fwidth)
        tax2 = ax2.twinx()
        setup_tsplot(ax2, xlim=(.1, 200))
        tax2.set_yscale('log')
        if args.verbose: 
            print('tb: %i, len(kmax): %i' % (tb.shape, len(kmax)))
        lns1 = ax2.plot(tb[1:], kmax, label=r'$k_{\textrm{max}}$', linewidth=2, color='black')
        ax2.set_ylabel(r'$k_{\textrm{max}}$')
        lns2 = tax2.plot(tb[1:], maxx, '--', color='black', linewidth=2, label='maxval(t)')
        ppo, ppcov = curve_fit(linear, np.log10(tb[11:]), np.log10(kmax[10:]))
        if args.verbose: 
            print(ppo)
        lns3 = tax2.plot(tb[1:], powerlaw(tb[1:], 10**(ppo[0]), ppo[1], x0=0.))
        lns = lns1 + lns2
        labs = [l.get_label() for l in lns]
        ax2.legend(lns, labs, loc='lower left', frameon=False)
        fig2.tight_layout(pad=0.3)
        savefig(fig2, args.ddir + '_peak.pdf')

def make_tsplot(fig, ax):
    ax.plot(tser.t, tser.brms, label=r'$B_{rms}$')
    ax.plot(tser.t, tser.urms, label=r'$u_{rms}$')
    ti = search_indx(tser.t, 10., eps=.05)
    if args.verbose: print(ti)
    if ti is not None:
        po1, pco1 = curve_fit(lambda x, a, b: powerlaw(x, a, b, x0=0), tser.t[ti:], tser.brms[ti:])
        po2, pco2 = curve_fit(lambda x, a, b: powerlaw(x, a, b, x0=0), tser.t[ti:], tser.urms[ti:])
        ax.plot(tser.t[10:], powerlaw(tser.t[10:],*po1))

def setup_tsplot(ax, xlim=None):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('time')
    if xlim is None:
        xlim = (1e-2, tser.t[-1])
    ax.set_xlim(xlim)

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

# read data
if args.spectra:
    tb, powerb = pc.read_power('power_mag.dat', datadir=args.ddir)
    if args.verbose: print('t.shape: ',tb.shape, ' last spectrum: ',tb[-1])
    p0 = (-15.56, 2.20, 4.26)

# Setup Figure and Axes
fwidth = 0.45
if args.peak:
    args.tsplot = True

if args.tsplot:
    fig, ax  = newfig(fwidth)
    tser = pc.read_ts(datadir=args.ddir, quiet=args.verbose)
    setup_tsplot(ax)
    make_tsplot(fig, ax)
 
if args.spectra:
    dim, krms = read_dim_krms(args.ddir)
    fig, ax  = newfig(fwidth)
    setup_specax(ax)
    make_specplot(fig, ax)
elif args.gridplot:
    fwidth = 0.95
    dim, krms = read_dim_krms(dirs[0])
    if len(dirs) == 4:
        rowcols = (2,2); r=2./3.
    elif len(dirs) == 6:
        rowcols = (2,3); r=0.45
    fig = plt.figure(figsize=figsize(.95, ratio=r))
    grid = Grid(fig, rect=[.08,.1,.9,.85], nrows_ncols=rowcols,
            axes_pad=0., label_mode='L',
            )
    #setup_gridax(ax)
    pb_ges, tb_ges = read_grid_data(dirs)
    strings = get_strings(dirs)
    for i,ax in enumerate(grid):
        plot_grid_spectra(ax, i)
        ax.title.set_visible(False)
        ax.text(1.5, 1e-4, strings[i] )
    fig.text(0.5, 0.01, 'k mode', ha='center')
    fig.text(0.01, 0.5, r'$E_{mag}$', va='center', rotation='vertical')

if not ( args.tsplot or args.spectra or args.gridplot):
    print('not making a plot, check settings...')
    sys.exit(1)

if not args.gridplot: 
    ax.legend(loc='lower center', ncol=2, frameon=False)
    fig.tight_layout(pad=.3)

if args.ddir.endswith('/'):
    args.ddir = args.ddir[:-1]
if args.spectra:
    plotname = '_mag_spec'
elif args.tsplot:
    plotname = '_tsplot'

if not args.gridplot:
    savefig(fig, args.ddir + plotname)
else:
    savefig(fig, 'grid_spec_' + args.ddir.split('_')[0])
print('SUCCESS')
