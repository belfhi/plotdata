import numpy as np
import pencil as pc
from os.path import join
import sys
ddir = sys.argv[1]
print(ddir)

Bx = np.load(join(ddir,'bx_shell.npy'))
By = np.load(join(ddir,'by_shell.npy'))
Bz = np.load(join(ddir,'bz_shell.npy'))

ux = np.load(join(ddir,'ux_shell.npy'))
uy = np.load(join(ddir,'uy_shell.npy'))
uz = np.load(join(ddir,'uz_shell.npy'))

jx = np.load(join(ddir,'jx_shell.npy'))
jy = np.load(join(ddir,'jy_shell.npy'))
jz = np.load(join(ddir,'jz_shell.npy'))

uxb1 = np.outer(uy, Bz) - np.outer(uz, By)
print('u x B _x : ',uxb1.shape)
uxb2 = np.outer(uz, Bx) - np.outer(ux, Bz)
uxb3 = np.outer(ux, By) - np.outer(uy, Bx)

T_kpq = np.tensordot(jx, uxb1,0) + np.tensordot(jy, uxb2,0) + np.tensordot(jz, uxb3,0) 
print(T_kpq.shape)
np.save('T_kpq.npy', T_kpq)

