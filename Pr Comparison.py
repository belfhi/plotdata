# coding: utf-8

import numpy as np
import matplotlib as mpl
mpl.use('pgf')

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

def newfig(width):
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename):
    plt.savefig('{}.pgf'.format(filename))
    plt.savefig('{}.pdf'.format(filename))

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


dirs = [s for s in listdir('.') if path.isdir(s) and not s.startswith('prandtl') ].sort()
print(dirs)

dim = pc.read_dim(datadir='pr=1/data')
krms = np.loadtxt('pr=1/data/power_krms.dat').flatten()[:dim.nxgrid//2]
print(krms.shape)
p0 = (-12.92758524, 1.94666781, 3.4643292)  #(-15.56, 2.20, 4.26)


# Simple plot
fig, ax  = newfig(0.45)
ax.set_xscale('log')
ax.set_xlim(.1,100)
ax.set_ylim(10,100)
for dd in dirs:   
    kmax = []
    if dd == 'pr=100':
        ddir = path.join(dd, 'old_data')
    elif dd == 'pr=1':
        ddir = path.join(dd, 'old2')
    else:
        ddir = path.join(dd, 'data')
    dim = pc.read_dim(datadir=ddir)
    if dd != 'pr=1':
        t, powerb = pc.read_power('power_mag.dat', datadir=ddir)
    else:
        t = np.loadtxt(path.join(dd, 'tb_resulting.dat'))
        powerb = np.loadtxt(path.join(dd, 'powerb_resulting.dat'))
        print(t.shape, powerb.shape)
    for p,pb in enumerate(powerb):
        xi = np.where( pb == pb.max())[0][0]; xi1 = xi - xi//2; xi2 = xi + xi//3
        po, pco = curve_fit(parabola, np.log10(krms[xi1:xi2]), np.log10(pb[xi1:xi2]), p0=p0)
        #if dd == 'pr=1':
        #    print(10**po[1])
        kmax.append(10**po[1])
    print(dd, len(kmax), t.shape)
    ax.plot(t, kmax, label='$%s$' %dd, linewidth=2.5)
ax.set_xlabel('time t')
ax.set_ylabel('$k_{max}$')
ax.set_yscale('log')
ax.set_yticks([10,20,50,100])
ax.set_yticklabels(['$10$','$20$','$50$','$100$'])
ax.legend(loc='lower center', ncol=4, prop={'size':20}, frameon=False)
fig.suptitle('Wavenumber of Maximum for different Prandtl Numbers', fontsize=24)

ax.set_xlabel('time')
ax.set_ylabel(r'$k_{\textrm{max}}$')

savefig('prandtl_scale_evolution')


# In[ ]:
#
#fig1, ax1 = plt.subplots(figsize=(12,9))
#tb, powerb = pc.read_power('power_mag.dat', datadir=dirs[0]+'/old2/')
#plottimes = []
#for t in [0, .1, 1, 10,  300]:
#    ti = search_indx(tb, t, eps= 0.01)
#    plottimes.append(ti)
#p0 = (-15.56, 2.20, 4.26)
#print(plottimes)
#kmax = []
#for p, pb in enumerate(powerb):
#    xi = np.where( pb == pb.max())[0][0]; xi1 = xi - xi//2; xi2 = xi + xi//3
#    #print(xi, xi1)
#    #print(xi)
#    po, pco = curve_fit(parabola, np.log10(krms[xi1:xi2]), np.log10(pb[xi1:xi2]), p0=p0)
#    if p in plottimes:
#        ax1.loglog(krms, pb, linewidth=2.5, label='t = %.0f' % tb[p])
#        ax1.plot(krms[xi1:xi2], 10**parabola(np.log10(krms[xi1:xi2]), *po), color='blue', linewidth=2.5)#, color=clrs[i])
#        print(10**po[1])
#    kmax.append(10**po[1])
#ax1.set_xlim(1,256)
#ax1.set_ylim(1e-11, 1e-3)
#ax1.legend(loc='lower center', ncol=2, frameon=False)
#
#
## In[ ]:
#
#fig = plt.figure(figsize=(12,7))
#ax = fig.add_subplot(111)
#ax.set_xscale('log')
#ax.set_xlim(.1,100)
#ax.set_ylim(10,100)
#for dd in dirs:   
#    kmax = []
#    if dd == 'pr=100':
#        ddir = path.join(dd, 'old_data')
#    elif dd == 'pr=1':
#        ddir = path.join(dd, 'old2')
#    else:
#        ddir = path.join(dd, 'data')
#    dim = pc.read_dim(datadir=ddir)
#    if dd != 'pr=1':
#        t, powerb = pc.read_power('power_mag.dat', datadir=ddir)
#    else:
#        t = np.loadtxt(path.join(dd, 'tb_resulting.dat'))
#        powerb = np.loadtxt(path.join(dd, 'powerb_resulting.dat'))
#        print(t.shape, powerb.shape)
#    for p,pb in enumerate(powerb):
#        xi = np.where( pb == pb.max())[0][0]; xi1 = xi - xi//2; xi2 = xi + xi//3
#        po, pco = curve_fit(parabola, np.log10(krms[xi1:xi2]), np.log10(pb[xi1:xi2]), p0=p0)
#        #if dd == 'pr=1':
#        #    print(10**po[1])
#        kmax.append(10**po[1])
#    print(dd, len(kmax), t.shape)
#    ax.plot(t, kmax, label='$%s$' %dd, linewidth=2.5)
#ax.set_xlabel('time t')
#ax.set_ylabel('$k_{max}$')
#ax.set_yscale('log')
#ax.set_yticks([10,20,50,100])
#ax.set_yticklabels(['$10$','$20$','$50$','$100$'])
#ax.legend(loc='lower center', ncol=4, prop={'size':20}, frameon=False)
#fig.suptitle('Wavenumber of Maximum for different Prandtl Numbers', fontsize=24)
#
#
## In[ ]:
#
#fig2 = plt.figure(figsize=(12,7))
#ax2 = fig2.add_subplot(111)
#ax2.set_xscale('log')
#ax2.set_yscale('log')
#ax2.set_xlim(.01,100)
#ax2.set_ylim(1e-2, .3)
#for dd in dirs:
#    ddir = path.join(dd, 'data')
#    dim = pc.read_dim(datadir=ddir)
#    tser = pc.read_ts(datadir=ddir)
#    ax2.plot(tser.t, tser.brms, label='$%s$' %dd, linewidth=2)
#ax2.set_xlabel('time t')
#ax2.set_ylabel('$B_{rms}$')
#ax2.legend(loc='lower left', ncol=2, prop={'size':20}, frameon=False)
#fig2.suptitle(r'Magnetic Field Strength for different Prandtl numbers $Pr=\eta/\nu$', fontsize=24)
#fig2.savefig('timeseries_comparison.pdf', bbox_inches='tight')
#
#
# In[1]:
