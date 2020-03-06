from params import *
import h5py
import pylab as pl


for i in range(t_eval.size):
    print('Time =', t_eval[i])

    h5f = h5py.File('fvm_data/%04d'%i + '.h5', 'r')
    sol = h5f['sol'][:]
    h5f.close()

    h5f = h5py.File('montecarlo_data/%04d'%i + '.h5', 'r')
    R   = h5f['R'][:]
    N   = h5f['N'][:]
    h5f.close()

    pl.plot(R * R0/1e-9, np.sum(sol[3:-3]) * dr * N, label = 'Particle Model')
    pl.plot(r * R0/1e-9, sol, 'k--', label = 'Distribution Model')
    
    pl.title('Time = %3.3f'%(t_eval[i] * t0/3600) + ' hours')
    pl.xlabel(r'$r$(in nm)')
    pl.ylabel(r'$f$')
    pl.xlim(0.4 * R0/1e-9, 2.5 * R0/1e-9)
    pl.ylim(0, 6)
    pl.legend(fontsize = 10)
    pl.savefig('comparison_images/%04d'%i + '.png', bbox_inches = 'tight')
    pl.clf()
