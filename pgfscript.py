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
from os import mkdir, path, listdir
from os.path import join, isdir, isfile

parser = argparse.ArgumentParser()
parser.add_argument('ddir', default='', type=str, help='datadir directory')
parser.add_argument('-s', '--spectra',  action='store_true', default=False, help="plot magnetic spectra")
parser.add_argument('-t', '--tsplot', action='store_true', default=False, help="plot timeseries too")
parser.add_argument('-g', '--gridplot', action='store_true', default=False, help='do the gridplot for multiple dirs')
parser.add_argument('-v', '--verbose', action='store_true', default=False, help='add verbosity')
parser.add_argument('-f', '--fit', action='store_true', default=False, help="don't fit a parabola")
parser.add_argument('-p', '--peak', action='store_true', default=False, help="make a plot of peak evolution")
parser.add_argument('-n', "--ncol", type=int, default=2, help="amount of columns in legend")
parser.add_argument('-d', "--demo", action='store_true', default=False, help="demonstrate the fitting method")
parser.add_argument('-k', "--kinetic", action='store_true', default=False, help="make kinetic spectra")
parser.add_argument('-l', "--light", action='store_true', default=False, help="Light color scheme")

args = parser.parse_args()
figdir = 'figures'

# need fit arg for peak:
if args.peak:
    args.fit = True

if args.ddir[-1] == '/':
    args.ddir = args.ddir[:-1]

def round_exp(num, prec=0):
    formatstr = '%.{0}e'.format(prec)
    newnum = float( formatstr % (num))
    return newnum

def check_dirname_grid(dd):
    if dd.split('_')[0] in ['visc', 'prandtl', 'slope','nonhel']:
        return True
    else:
        return False

if args.verbose:
    for arg in vars(args):
        print('{:<9s}: {:<12s}'.format(arg, str(getattr(args, arg))))

if args.gridplot and not check_dirname_grid(args.ddir):
    print('gridplot only works for nu comparison and prandtl numbers')
    print('usage:')
    print('   python pgfscript.py -g visc_')
    print(' or:')
    print('   python pgfscript.py -g prandtl')
    print(' or:')
    print('   python pgfscript.py -g slope')
    sys.exit(1)
elif args.gridplot:
    dirs = [d for d in listdir('.') if d.startswith(args.ddir[:4]) and isdir(d) and not isfile(join(d,'NOPLOT'))]
    try:
        dirs = sorted(dirs, key=lambda a: float(a.split('_')[-1]))
    except ValueError:
        dirs = sorted(dirs, key=lambda a: float(a.split('_')[-1][1:]))
    if args.verbose: print(dirs)

def check_figdir(directory):
    if not path.exists(directory):
        mkdir(directory)
        if path.isdir(directory):
            return
    else:
        return
    
def read_power_data(ddir):
    if args.ddir[:4] == 'zeus':
        tb, powb = read_zeus_data(args.ddir)
    else:
        try:
            powb = np.loadtxt(path.join(ddir,'powertot.dat'))
            tb = np.loadtxt(path.join(ddir,'ttot.dat'))
            if args.verbose: print('reading powertot.dat File')
        except FileNotFoundError:
            tb, powb = pc.read_power('power_mag.dat', datadir=ddir)
    if args.verbose: print('t.shape: %i, last spectrum: %.f' %(tb.shape[0], tb[-1]))

    return tb, powb
    

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
    return r'{0}\times 10^{{{1}}}'.format(s1, s2)

def setup_specax(ax):
    ax.set_yscale('log')
    ax.set_xscale('log')
    i = 1 if args.ddir.startswith('delta') else 0
    ax.set_xlim(1,dim.nxgrid//2)
    ax.set_xlabel('k mode')
    if args.kinetic:
        ax.set_ylabel(r'$E(k,t)$')
        ax.text(.05,.95, 'solid: magnetic\ndashed: kinectic', transform=ax.transAxes, 
                horizontalalignment='left', verticalalignment='top')
    else:  
        ax.set_ylabel(r'$E_B(k,t)$')
    
    #ylims = (round_exp(powerk[i+1,1]), round_exp( 1.5*powerk[i,:].max()))
    ylims = (round_exp(powerb[i+1,1]), round_exp( 1.5*powerb[i,:].max()))
    ax.set_ylim(ylims)
    if args.verbose:
        print('ylim : %.1e %.1e ' % (ax.get_ylim()))

def get_plottimes(tb):
    ptimes = []
    if args.ddir.startswith('delta'):
        pptimes = [0,.1, 1, 10,200]
    elif args.ddir.startswith('nonhel') and tb[-1] > 300:
        pptimes = [1,5,50,300]
    elif args.gridplot and args.ddir.startswith('prandtl'):
        pptimes = [0,1,10,50, 300]
    elif args.gridplot and args.ddir.startswith('visc'):
        pptimes = [0,1,10,50, 100]
    elif args.gridplot and args.ddir.startswith('slope'):
        pptimes = [0,1,10,100]
    elif args.gridplot and args.ddir.startswith('nonhel'):
        pptimes = [0,1,10,100]
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

def plot_spectra(ax, ddir):
    tb, powerb = read_power_data(ddir)
    plottimes = get_plottimes(tb)
    clrindx = iter(np.linspace(0,1,len(plottimes)))
    for p, pb in enumerate(powerb):
        if p in plottimes:
            clr = cscheme(next(clrindx))
            if args.verbose: print('clr: %.1f, %.1f, %.1f, %.1f' % clr)
            ax.loglog(krms, pb, linewidth=1.5, color=clr)
    ax.set_xlim(1,dim.nxgrid//2)
    ymax = round_exp(3*powerb[0,:].max())
    #global ymax0
    #if ymax > ymax0:
    ax.set_ylim(powerb[0,1], ymax)
    ymax0 = ymax
    return ax

def make_specplot(fig, ax, powerarr, kin=False):
    p0 = (-15.56, 2.20, 4.26)
    plottimes = get_plottimes(tb)
    clrindx = iter(np.linspace(0,1,len(plottimes)))
    kmax = []
    maxx = []
    for p, pb in enumerate(powerarr):
        if not (p == 0 and args.ddir.startswith('delta')):
            xi = np.where( pb == pb.max())[0][0]; xi1 = xi - xi//3; xi2 = xi + xi//2
            #if args.verbose: print('xi: %i, xi1: %i, xi2: %i' %(xi, xi1, xi2))
            if args.fit:
                po, pco = curve_fit(parabola, np.log10(krms[xi1:xi2]), np.log10(pb[xi1:xi2]), p0=p0)
                kmax.append(10**(po[1]))
                maxx.append(max(pb))
        else:
            ax.plot(krms[1:], pb[1:], color=cscheme(next(clrindx)), linewidth=1.5)
            continue
        #print(xi, xi1, xi2)
        if p in plottimes:
            if tb[p] < 1.0 and tb[p] > 0:
                s = 't = %.1f' 
            else:
                s = 't = %.0f'
            #s = 't = %.0f'
            curr_clr = cscheme(next(clrindx))
            if not kin:
                ax.plot(krms[1:], pb[1:], label=s % tb[p], color=curr_clr, linewidth=1.5)
            else:
                ax.plot(krms[1:], pb[1:], '--', color=curr_clr, linewidth=1.5)
            if args.verbose: print('t= %.1f, color: %.2f, %.2f, %.2f, %.2f' % (tb[p], *curr_clr))
            if args.demo:
                ax.plot(krms[xi1], pb[xi1], "o", color='blue', linewidth=1.5)
                ax.plot(krms[xi2], pb[xi2], "o", color='blue', linewidth=1.5)
                xvals = np.arange(krms[xi1],krms[xi2],.1)
                ax.plot(xvals, 10**parabola(np.log10(xvals), *po), color='blue')
    if args.peak:
        fig2, ax2 = newfig(fwidth)
        tax2 = ax2.twinx()
        setup_tsplot(ax2, xlim=(tb[1], 300))
        print(ax2.get_xlim())
        tax2.set_yscale('log')
        if args.verbose: 
            print('tb: %i, len(kmax): %i' % (tb.shape[0], len(kmax)))
        i = 1 if args.ddir.startswith('delta') else 0
        ppo, ppcov = curve_fit(linear, np.log10(tb[i+10:]), np.log10(kmax[10:]))
        ppo1, ppcov1 = curve_fit(linear, np.log10(tb[i+1:]), np.log10(maxx[1:]))
        lns4 = tax2.plot(tb[1:], powerlaw(tb[1:], 10**(ppo1[0]), ppo1[1], x0=0.), color=clr1,
                         label=r'fit: $E(k_{max},t)\sim t^{%.0f}$' % ppo1[1], linewidth=2)
        lns3 = ax2.plot(tb[1:], powerlaw(tb[1:], 10**(ppo[0]), ppo[1], x0=0.), linewidth=2,
                        label=r'fit: $k_{max}\sim t^{%.1f}$' % ppo[1], color=clr2)
        lns1 = ax2.plot(tb[i:], kmax, label=r'$k_{\textrm{max}}$', linewidth=2, color='black')
        lns2 = tax2.plot(tb[i:], maxx, '--', color='black', linewidth=2, label='Spectrum Max value')
        ax2.set_ylabel(r'$k_{\textrm{max}}$')
        #tax2.set_ylabel('max(E(k,t))')
        if args.verbose: 
            print(ppo, ppo1)
        lns = lns1 + lns3 + lns2 + lns4
        labs = [l.get_label() for l in lns]
        ax2.legend(lns, labs, loc='lower left', ncol=1, frameon=False)
        fig2.tight_layout(pad=0.3)
        savefig(fig2, path.join(figdir,'%s_peak' % args.ddir) )

def make_tsplot(fig, ax):
    ax.plot(tser.t, tser.brms**2, label=r'$E_{\textrm{mag}}$', linewidth=1.5, color=clr1)
    ax.plot(tser.t, tser.urms**2, label=r'$E_{\textrm{kin}}$', linewidth=1.5, color=clr2)
    ti = search_indx(tser.t, 10., eps=.05)
    ti2 = search_indx(tser.t, 200., eps=.05)
    if args.verbose: print('index from which fitting begins: ',ti)
    if ti is not None:
        po1, pco1 = curve_fit(lambda x, a, b: powerlaw(x, a, b, x0=0), tser.t[ti:], tser.brms[ti:]**2)
        po2, pco2 = curve_fit(lambda x, a, b: powerlaw(x, a, b, x0=0), tser.t[ti:], tser.urms[ti:]**2)
        ax.plot(tser.t[ti:ti2], powerlaw(tser.t[ti:ti2],1.5*po1[0],po1[1]), linewidth=1.5, 
                label=r'fit: $B\sim t^{%.1f}$' % po1[1], color=clr3)
        if args.verbose:
            print('Magnetic Energy: E={0:.1e}*t^{1:.1f}'.format(*po1))
            print('Kinetic Energy: E={0:.1e}*t^{1:.1f}'.format(*po2))

def setup_tsplot(ax, xlim=None):
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('time')
    ax.set_ylabel(r'$E_{B},\ E_{K}$')
    if xlim is None:
        xlim = (1e-2, tser.t[-1])
    ax.set_xlim(xlim)

def read_zeus_data(ddir):
    files = sorted([s for s in listdir(ddir) if s.startswith('magspec_')], 
            key=lambda s: int(s.split('_')[1][:-4]))
    if args.verbose:
        print(files, type(files))
    tb = np.zeros(len(files))
    powerb = np.zeros((len(files),256))
    for f, fi in enumerate(files):
        with open(join(ddir, fi)) as infile:
                tb[f] = int(float(infile.readline())) 
                print('read spec at: t=', tb[f])
                powerb[f,:] = [float(x) for x in infile.readlines()]
    return tb, powerb

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

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5
right = left + width
top = bottom + height

# read kinetic data
if args.kinetic:
    tk, powerk = pc.read_power('power_kin.dat', datadir=args.ddir)
    if args.verbose:
        print('tk shape: ', tk.shape, 'last spectrum: ', tk[-1])

# Setup Figure and Axes
fwidth = 0.45
#if args.peak:
#    args.tsplot = True

if args.tsplot:
    fig, ax  = newfig(fwidth)
    if args.verbose:
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        print('ts figure width: %.1f, height: %.1f' % (width, height))
    tser = pc.read_ts(datadir=args.ddir, quiet=not args.verbose)
    setup_tsplot(ax)
    make_tsplot(fig, ax)
 
if args.spectra:
    tb, powerb = read_power_data(args.ddir)
    dim, krms = read_dim_krms(args.ddir)
    fig, ax  = newfig(fwidth)
    if args.verbose:
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        print('spec figure width: %.1f, height: %.1f' % (width, height))

    setup_specax(ax)
    make_specplot(fig, ax, powerb, kin=False)
    if args.kinetic:
        make_specplot(fig, ax, powerk, kin=True)
elif args.gridplot:
    fwidth = 0.95
    ymax0 = 0.
    dim, krms = read_dim_krms(dirs[0])
    if len(dirs) == 4:
        rowcols = (2,2); r=2/3
    elif len(dirs) == 6:
        rowcols = (2,3); r=0.45
    elif len(dirs) == 8:
        rowcols = (2,4); r=0.33
    else:
        rowcols = (len(dirs)//2, 2); r=1.33
        print('warning: guessing row / column combination')
    fig = plt.figure(figsize=figsize(.95, ratio=r))
    grid = Grid(fig, rect=[.08,.1,.9,.85], nrows_ncols=rowcols,
                axes_pad=0., label_mode='L' )
    #setup_gridax(ax)
    #pb_ges, tb_ges = read_grid_data(dirs)
    strings = [get_string(d) for d in dirs]
    for i,ax in enumerate(grid):
        plot_spectra(ax, dirs[i])
        ax.title.set_visible(False)
        # axes coordinates are 0,0 is bottom left and 1,1 is upper right
        ax.text(.05, .95, strings[i], horizontalalignment='left',
                verticalalignment='top',transform=ax.transAxes)
        #ax.text(1.5, y1, strings[i] )
    fig.text(0.5, 0.01, 'k mode', ha='center')
    fig.text(0.01, 0.5, r'$E_{\textrm{mag}}$', va='center', rotation='vertical')

if not ( args.tsplot or args.spectra or args.gridplot):
    print('not making a plot, check settings...')
    sys.exit(1)

if not args.gridplot: 
    if args.tsplot:
        ax.legend(loc='lower left', ncol=1, frameon=False)
        if args.verbose: print('legend low left')
    else:
        ax.legend(loc='lower center', ncol=args.ncol, frameon=False)
    fig.tight_layout(pad=.3)

if args.ddir.endswith('/'):
    args.ddir = args.ddir[:-1]
if args.spectra and not args.kinetic:
    plotname = '_mag_spec'
elif args.spectra and args.kinetic:
    plotname = '_kin_spec'
elif args.tsplot:
    plotname = '_tsplot'

if not args.gridplot:
    filename = args.ddir + plotname
else:
    filename = 'grid_spec_' + args.ddir.split('_')[0]
if args.verbose:
    print(filename)
savefig(fig, path.join(figdir, filename))
print('SUCCESS')
