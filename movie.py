from IPython.display import HTML, Latex, Javascript
from matplotlib import animation
import matplotlib as mpl
mplpars = {'mathtext.fontset':'stix', 'font.family': 'STIXGeneral', 'font.size': 24,
           'lines.linewidth': 2, 'figure.dpi': 100, 'figure.figsize': (10/3,2.5),
           'legend.fontsize': 22.0, 'legend.frameon': False}
mpl.rcParams.update(mplpars)
import pencil as pc
import numpy as np
import matplotlib.pyplot as plt
from helperfuncs import search_indx
from os.path import join

ddir = 'reionization/k-0.9/'
tb, powerb = pc.read_power('power_mag.dat', datadir=ddir)
dim = pc.read_dim(datadir=ddir)
krms = np.loadtxt(join(ddir, 'power_krms.dat')).flatten()[:dim.nxgrid//2]
tser = pc.read_ts(datadir=ddir)
print(tb[-1], tb.shape, powerb.shape)

fig = plt.figure(figsize=(9,7))
ax = plt.axes(xlim=(1,768), ylim=(1e-7, 1e-2))#xlim=(1,512), ylim=(1e-7, 1e-2))
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('$k$')# mode')
ax.set_ylabel(r'$E_B(k)$')
line, = ax.plot([], [], lw=2)
#line2, = ax.plot([], [], lw=2)
txt = ax.text(7e2,9e-3, 'z = ', va='top', ha='right')
fig.tight_layout()
# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    #line2.set_data([], [])
    return line,#line2

# animation function.  This is called sequentially
def animate(i):
    ts = search_indx(tser.t, tb[i], eps=0.01)
    powi = powerb[i,1:]
    #fiti = logcutoff(krms[1:], *opts[i,:])
    txt.set_text('z = %.0f' %tser.zredshift[ts])
    #line.set_data(np.log10(krms[1:]), fiti)
    line.set_data(krms[1:], powi)
    #line.set_data(x,y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=130, interval=100, blit=True, repeat=False)

#plt.close()
AVwrite = animation.MencoderWriter()
fname = 'anim_%s.mp4'# % ddir[:-1]
anim.save(fname , writer=AVwrite)
# call our new function to display the animation

