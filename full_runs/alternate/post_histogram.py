from params import *
import h5py
import pylab as pl

maximum = 0
for i in range(t_eval.size):
    h5f = h5py.File('montecarlo_data/%04d'%i + '.h5', 'r')
    R   = h5f['R'][:]
    N   = h5f['N'][:]
    h5f.close()

    if(N.max() > maximum):
        maximum = N.max()

for i in range(t_eval[:251].size):
    print('Time =', t_eval[i])

    h5f = h5py.File('montecarlo_data/%04d'%i + '.h5', 'r')
    R   = h5f['R'][:]
    N   = h5f['N'][:]
    h5f.close()

    R = R * R0/1e-9
    width = 0.9 * (R[1] - R[0])

    pl.bar(R, N, align='center', width=width)
    pl.title('Time = %3.3f'%(t_eval[i] * t0/3600) + ' hours (growth)')	
    pl.xlabel(r'$r$ (in nm)')
    pl.ylabel(r'$N$')
    pl.xlim(0.4 * R0/1e-9, 2.5 * R0/1e-9)
    pl.ylim(0, maximum + 10)
    pl.savefig('histogram_images1/%04d'%i + '.png', bbox_inches = 'tight', dpi = 300)
    pl.clf()

for i in range(251, t_eval.size):
    print('Time =', t_eval[i])

    h5f = h5py.File('montecarlo_data/%04d'%i + '.h5', 'r')
    R   = h5f['R'][:]
    N   = h5f['N'][:]
    h5f.close()

    R = R * R0/1e-9
    width = 0.9 * (R[1] - R[0])

    pl.bar(R, N, align='center', width=width)
    pl.title('Time = %3.3f'%(t_eval[i] * t0/3600) + ' hours (ripening)')
    pl.xlabel(r'$r$ (in nm)')
    pl.ylabel(r'$N$')
    pl.xlim(0.4 * R0/1e-9, 2.5 * R0/1e-9)
    pl.ylim(0, maximum + 10)
    pl.savefig('histogram_images2/%04d'%(i-251) + '.png', bbox_inches = 'tight', dpi = 300)
    pl.clf()
