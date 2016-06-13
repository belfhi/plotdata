#!/usr/bin/env python3
import numpy as np
import matplotlib as mpl
from os import path, listdir
mpl.use('pgf')

def figsize(scale, ratio=None):
    fig_width_pt = 523.5307                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    if not ratio:
        golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
        fig_height = fig_width*golden_mean              # height in inches
    else:
        fig_height = fig_width*ratio
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
    "text.fontsize": 10,
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

# I make my own newfig and savefig functions
def newfig(width):
    plt.clf()
    fig = plt.figure(figsize=figsize(width))
    ax = fig.add_subplot(111)
    return fig, ax

def savefig(filename):
    plt.savefig('{}.pgf'.format(filename))
    plt.savefig('{}.pdf'.format(filename))


# Simple plot
fig, ax1  = newfig(0.6)


specfiles = sorted(filter(lambda s: s.startswith('magspec_'), listdir('.')))
print(specfiles)

for magsp in specfiles:
    with open(magsp) as magf:
        t = float(magf.readline())
        #print(t)
        if t == 0.0:
            t += 0.1
        spec = [float(x) for x in magf.readlines()]
    if round(t,1) in [0.1, 1, 5, 10, 100]:
        s = 't = %.1f'  if t <= 0.1 else 't = %.0f'
        ax1.plot( spec, label=s%t, linewidth=1.5)    spec = []
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('k mode')
ax1.set_ylabel('$E_{mag}(k)$')
ax1.set_xlim(1,256)
ax1.set_ylim(1e-9, 1e-4)
ax1.plot(range(20,100), [1e-3*x**-2 for x in range(20,100)], color='black', linewidth=2)
ax1.text(30,1e-7, '$\propto k^{-2}$')
ax1.legend(loc='center left', bbox_to_anchor=(1.,.5), prop={'size':16}, frameon=False)
fig.tight_layout(pad=0.3)
#fig.suptitle('Zeus-MP2 Run, $512^3$', fontsize=24)
savefig('zeus_mp2_pencil_init')
#ff = h5py.File('Combined/hdfaa.010')
#mag1 = ff.get('i_mag_field')
#mag2 = ff.get('j_mag_field')
#mag3 = ff.get('k_mag_field')
#ti = ff.get('time')
#t0 = ti[0]
#print(t0)
##vel1 = ff.get('i_velocity')
##dens = ff.get('Density')
##dd = dens[...]
##v1 = vel1[...]
##print(mm1.max(), np.where(mm1 == mm1.max()))#, dd[np.where(mm1 == mm1.max())])
##print(dd.max(), np.where(dd == dd.max()), v1[np.where(dd == dd.max())])
###print(v1.max(), np.where(v1 == v1.max()), dd[np.where(v1 == v1.max())])
#
#absmag = np.sqrt( mag1[255,:,:]**2 + mag2[255,:,:]**2 + mag3[255,:,:]**2)
#
#fig1 = plt.figure(figsize=(9,9))
#ax1 = fig1.add_subplot(111)
#im1 = ax1.imshow(absmag[:,:], cmap=plt.cm.hot, extent=[-np.pi, np.pi, -np.pi, np.pi])
#pis = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
#pistr = ['$-\pi$', '$-\pi/2$', '$0$', '$\pi/2$','$\pi$']
#ax1.set_xticks(pis)
#ax1.set_xticklabels(pistr)
#ax1.set_yticks(pis)
#ax1.set_yticklabels(pistr)
#finame = 'zeus_init_slice_t%.0f.png' %t0
#print(finame)
#fig1.savefig(finame, bbox_inches='tight')
