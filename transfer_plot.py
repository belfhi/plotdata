#!/usr/bin/env python3
import matplotlib as mpl
mplpars = {'mathtext.fontset':'stix', 'font.family': 'STIXGeneral', 'font.size': 20,
           'lines.linewidth': 2, 'figure.dpi': 72, 'figure.figsize': (10,8),
           'legend.fontsize': 20.0, 'legend.frameon': False,}
mpl.rcParams.update(mplpars)
import matplotlib.pyplot as plt
import numpy as np
import pencil as pc
import sys
from scipy.ndimage.filters import gaussian_filter
from matplotlib.ticker import MultipleLocator, Formatter
import argparse
from os.path import join, exists

parser = argparse.ArgumentParser()
parser.add_argument('ddir', type=str, default='nonhel_k80', help='select datadir to read files from')
parser.add_argument('-s', '--slice',  action='store', type=int, default=0, choices=[0,1,2], help="select if you want to plot T_k(p,p), T_q(k,p) or T_p(q,k)")
parser.add_argument('-i', '--index', action='store', default=10, type=int, help="plot timeseries too") 
args = parser.parse_args()

#read T_kpy.npy file
try:
    T_kpq = np.load(join(args.ddir, 'T_kpq.npy'))
    print('T_kpq shape: ', T_kpq.shape)
except FileNotFoundError:
    print('no T_kpy.npy file found. Aborting.')
    sys.exit(1)

fig1, ax1 = plt.subplots( figsize=(5,5))
ti = MultipleLocator(base=1.0)
i = args.index
j = args.slice
ex = 512
sigma = 2
k0=80

xlbl = ['$p/p_0$','$q/q_0$', '$k/k_0$']
ylbl = ['$q/q_0$','$k/k_0$', '$p/p_0$']
ax1.text(0.9,0.05, '{0} = {1!s}'.format(xlbl[j-1], i/k0), 
         transform=ax1.transAxes, ha='right')
Tk = T_kpq[i,:,:]
Tp = T_kpq[:,i,:]
Tq = T_kpq[:,:,i]

Tk_smooth =0.5*( Tk[::2,::2] + Tk[1::2,1::2])
Tp_smooth =0.5*( Tp[::2,::2] + Tp[1::2,1::2])
Tq_smooth =0.5*( Tq[::2,::2] + Tq[1::2,1::2])

T_list = [Tk, Tp, Tq]
T_smooth = [Tk_smooth, Tp_smooth, Tq_smooth]
gauss_image = gaussian_filter(T_list[j], sigma)
#ax2.imshow(np.real(T_kpq[256,:,:]), cmap=plt.cm.hot)
im1 = ax1.imshow(gauss_image[:ex,:ex], extent=[0,ex/k0,ex/k0,0])
#im2 = ax2.imshow(T_list[j][:ex,:ex], extent=[0,ex/k0,ex/k0,0])
#divider = make_axes_locatable(ax2)
#ax2.set_xlim(0,100/k0)
#ax2.set_ylim(0,100/k0)
#cax = divider.append_axes("right", size="5%", pad=0.05)
ax1.set_xlabel(xlbl[j])
ax1.set_ylabel(ylbl[j])
ax1.yaxis.set_major_locator(ti)
ax1.xaxis.set_major_locator(ti)
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top')
#plt.colorbar(im, cax=cax)
fform = 'png'
fname = 'T_kpq_plot{0}.{1}'.format(str(j),fform)
print(fname)
fig1.tight_layout()
fig1.savefig(fname)
