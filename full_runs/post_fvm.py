import numpy as np
import h5py
import pylab as pl
pl.style.use('latexplot')

r_min = 0.5
r_max = 5
N_r   = 256
N_g   = 3
dr    = (r_max - r_min) / N_r

r  = r_min + (0.5 + np.arange(-N_g, N_r + N_g)) * dr

h5f = h5py.File('fvm_data/0001.h5', 'r')
sol = h5f['sol'][:]
h5f.close()

pl.plot(r, sol[0])
pl.title('Time = 000' + r'$\tau$')
pl.xlabel(r'$r$')
pl.ylabel(r'$f$')
pl.ylim(0, 5)
pl.savefig('fvm_images/0000.png', bbox_inches = 'tight')
pl.clf()

t_final = 240
for T in range(1, t_final):
    if(T%10 == 0):
        print('Time =', T)
        h5f = h5py.File('fvm_data/%04d'%T + '.h5', 'r')
        sol = h5f['sol'][:]
        h5f.close()

        pl.plot(r, sol[-1])
        pl.title('Time = %03d'%T + r'$\tau$')
        pl.xlabel(r'$r$')
        pl.ylabel(r'$f$')
        pl.ylim(0, 5)
        pl.savefig('fvm_images/%04d'%T + '.png', bbox_inches = 'tight')
        pl.clf()
