import matplotlib as mpl
mpl.rcParams.update({'font.size': 20, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})
mpl.rcParams.update({'lines.linewidth': 2, 'figure.dpi': 100, 'legend.fontsize': 20.0, 'legend.frameon': False})
import matplotlib.pyplot as plt
import pencil as pc
import numpy as np
from mpipencil import read_slices, read_vector_slices
from os import path, mkdir
#
xyticks = [-np.pi,-np.pi/2,0,np.pi/2,np.pi]
texticks = ['$-\pi$','$-\pi/2$','$0$', '$\pi/2$','$\pi$']
#
figdir = 'figures'
try:
    mkdir(figdir)
except FileExistsError:
    pass
def set_piticks(axes):
    axes.set_yticks(xyticks)
    axes.set_yticklabels(texticks)
    axes.set_xticks(xyticks)
    axes.set_xticklabels(texticks)
    #return axes
#
ddir = 'helical'
bb_slices = read_vector_slices(datadir=ddir, field='bb', absval=True)
print(bb_slices.nslice, bb_slices.t[0], bb_slices.t[-1])
#
fig, ax = plt.subplots(figsize=(4,4))
txt = ax.text(-np.pi+0.3,np.pi-.3, '', ha='left', va='top',
        bbox={'boxstyle':'square, pad=0.3', 'fc':'white', 'alpha':0.5})
x = [-np.pi, np.pi, -np.pi, np.pi]
set_piticks(ax)
for i in range(4,bb_slices.nslice,5):
    ti = bb_slices.t[i]
    if i == 4:
        ti = ti-0.5
    print(i, ti )
    txt.set_text('t=%.0f' %ti)
    ax.imshow(bb_slices.bb[i,:,:], extent=x, cmap=plt.cm.hot)
    fig.savefig(path.join(figdir,'bb_slice_%03i.pdf') %(i-4), bbox_inches='tight')
    fig.savefig(path.join(figdir,'bb_slice_%03i.png') %(i-4), bbox_inches='tight')
