import numpy as np
import h5py
import pylab as pl
pl.style.use('latexplot')

r_min = 0.5
r_max = 3
N_r   = 256
N_g   = 3
dr    = (r_max - r_min) / N_r

r  = r_min + (0.5 + np.arange(-N_g, N_r + N_g)) * dr

h5f = h5py.File('fvm_data/0001.h5', 'r')
sol = h5f['sol'][:]
h5f.close()

h5f = h5py.File('montecarlo_data/0000.h5', 'r')
R = h5f['R'][:]
N = h5f['N'][:]
h5f.close()

pl.plot(R, N)
pl.plot(r, sol[0], 'k--')
pl.title('Time = 000' + r'$\tau$')
pl.xlabel(r'$r$')
pl.ylabel(r'$f$')
pl.xlim(0.4, 2)
pl.ylim(0, 6)
pl.savefig('comparison_images/0000.png', bbox_inches = 'tight')
pl.clf()

t_final = 1000
for T in np.arange(1, t_final):
    print('Time =', T)

    h5f = h5py.File('fvm_data/%04d'%T + '.h5', 'r')
    sol = h5f['sol'][:]
    h5f.close()

    h5f = h5py.File('montecarlo_data/%04d'%T + '.h5', 'r')
    R   = h5f['R'][:]
    N   = h5f['N'][:]
    h5f.close()

    pl.plot(R, N)
    pl.plot(r, sol[-1], 'k--')
    pl.title('Time = %03d'%T + r'$\tau$')
    pl.xlabel(r'$r$')
    pl.ylabel(r'$f$')
    pl.xlim(0.4, 2)
    pl.ylim(0, 6)
    pl.savefig('comparison_images/%04d'%T + '.png', bbox_inches = 'tight')
    pl.clf()
