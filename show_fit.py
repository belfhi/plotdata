# coding: utf-8

import numpy as np
import matplotlib as mpl
import sys
from os import path, listdir
mpl.use('pgf')

try:  
    dstart = sys.argv[1]
    tfit = int(sys.argv[2])
except IndexError:
    print('must give directory and fitting index')
    sys.exit(1)

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

def round_exp(num, prec=0):
    formatstr = '%.{0}e'.format(prec)
    newnum = float( formatstr % (num))
    return newnum

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
    return r'{0}{{\times}} 10^{{{1}}}'.format(s1, s2)

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
from scipy.optimize import curve_fit

try:
    dim = pc.read_dim(datadir=dstart)
    krms = np.loadtxt(path.join(dstart,'power_krms.dat')).flatten()[:dim.nxgrid//2]
    print('krms shape: ', krms.shape)
    p0 = (-12.92758524, 1.94666781, 3.4643292)  #(-15.56, 2.20, 4.26)
except FileNotFoundError:
    sys.exit(1)

#read in simulation parameters 
dim = pc.read_dim(datadir=dstart)
pars = pc.read_param(datadir=dstart, quiet=True, param2=True)
xi = dim.nxgrid//2
t, powerb = pc.read_power('power_mag.dat', datadir=dstart)
print('t: ',t.shape, 'powerb: ',powerb.shape)
print('plotting at t = %g' % t[tfit])

# Simple plot setup
plotted = False
fig, ax  = newfig(0.45, ratio=0.75)
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xlim(1,dim.nxgrid//2)
ax.set_xlim(10,50)
ax.set_ylim(round_exp(powerb[tfit,10]),1.5*max(powerb[tfit,:]))
ax.set_xticks([10,20,50])
ax.set_xticklabels([str(s) for s in ax.get_xticks()])
#ax.set_ylim(1e-7, 1e-3)
clr = (0.9, 0.4, 0.4, 1.0)
for p,pb in enumerate(powerb):
    try:
        if t[p] < pars.tforce_stop:
            continue
    except AttributeError:
        pass
    xi = np.where( pb == pb.max())[0][0]; xi1 = xi - xi//3; xi2 = xi + xi//2
    if xi2 > dim.nxgrid//2:
        continue
    #print(xi, xi1, xi2)
    xvals = np.arange(krms[xi1], krms[xi2], 0.1)
    try:
        po, pco = curve_fit(parabola, np.log10(krms[xi1:xi2]), np.log10(pb[xi1:xi2]), p0=p0)
    except RuntimeError:
        print(xi, xi1, xi2)
        continue
    if p == tfit:
        plotted=True
        ax.plot(xvals, 10**(parabola(np.log10(xvals), *po)), color=clr, linewidth=2, label='fitted parabola')
        ax.plot(krms[1:], pb[1:], color='black', label='simulation data')
        ax.plot(krms[xi1], pb[xi1], 'o', color=clr)
        ax.plot(krms[xi2], pb[xi2], 'o', color=clr)
    p0 = po
#print(len(kmax), t.shape)
#ax.plot(t[1:], kmax[1:], label=s  , linewidth=1.5)
#ax.set_yticks([10,20,50,100])
#ax.set_yticklabels(['$10$','$20$','$50$','$100$'])
ax.legend(loc='lower center', ncol=1, frameon=False)
#fig.suptitle('Wavenumber of Maximum ')

ax.set_xlabel('k mode')
ax.set_ylabel(r'$E_{\textrm{mag}}$')
fig.tight_layout(pad=0.3)
fname = 'fitting_example_%s' % dstart
filename = path.join('figures',fname)
print(filename)
savefig(filename)
