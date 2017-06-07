#!/bin/env python3

from mpi4py import MPI
import numpy as np
import pencil as pc
import numpy.fft as F
from mpipencil import read_var
import sys
from os.path import join, exists

def integrate_shell(fft_arr, d):
    shell = np.zeros(shape=(d.nx//2))
    for iz in range(d.nz):
        for iy in range(d.ny):
            for ix in range(d.nx):
                k = int(round(np.sqrt(kx[ix]**2  + kx[ix]**2 + kz[iz]**2)))
                if (k>= 0 and k<=nk-1):
                    shell[k] += np.sqrt(np.real(fft_arr[ix,iy,iz])**2 + np.imag(fft_arr[ix,iy,iz])**2)
    return shell

folder = 'field_arr_files'
ddir = 'data'
comm = MPI.COMM_WORLD
print('initialized MPI')
size = comm.Get_size()
rank = comm.Get_rank()
if rank == 0:
    print('size = ',size)
files = ['bx.npy', 'by.npy', 'bz.npy', 'jx.npy', 'jy.npy', 'jz.npy', 'ux.npy', 'uy.npy', 'uz.npy']
lfile = files[rank]
lfield = lfile[:2]
print(lfield)
if len(files) != size:
    comm.Abort()

fname = join(folder, lfile)
print(fname)
try:
    arr = np.load(fname)
except FileNotFoundError:
    print('no such file')
    comm.Abort()

dim = pc.read_dim(datadir=ddir)
nk = dim.nx//2
if rank == 0:
    print('nk = ',nk)
kx = np.roll(np.arange(-nk,nk), nk)
ky = np.roll(np.arange(-nk,nk), nk)
kz = np.roll(np.arange(-nk,nk), nk)
arr_fft = F.fftn(arr, norm='ortho')
print(arr_fft.shape, arr_fft.max())
shell_int = integrate_shell(arr_fft, dim)
fname_shell = lfield + '_shell.npy'
np.save(fname_shell, shell_int)
