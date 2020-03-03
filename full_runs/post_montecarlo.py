import numpy as np
import h5py
import pylab as pl
pl.style.use('latexplot')

t_final = 200
for T in np.arange(0, t_final):
    print('Time =', T)
    h5f = h5py.File('montecarlo_data/%04d'%T + '.h5', 'r')
    R = h5f['R'][:]
    N = h5f['N'][:]
    h5f.close()

    pl.plot(R, N)
    pl.title('Time = %03d'%T + r'$\tau$')
    pl.xlabel(r'$r$')
    pl.ylabel(r'$f$')
    pl.xlim(0, 5)
    pl.ylim(0, 5)
    pl.savefig('montecarlo_images/%04d'%T + '.png', bbox_inches = 'tight')
    pl.clf()
