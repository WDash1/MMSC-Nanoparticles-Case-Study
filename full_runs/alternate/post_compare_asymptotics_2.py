from params import *
import h5py
import pylab as pl


def analytical_1(t):
    A = cs * np.exp(lcap / (R0*r_min))
    s = np.sqrt((r + Da)**2 - 2*(cinf0*t - A*t) / delta_C) - Da
    f = np.sqrt(0.5/np.pi) / sigma * np.exp(-1/2 * ((s - 1)/sigma)**2)
    return f



def analytical_2(t):
    f_prime = cs*np.exp(lcap/(R0*r))/(delta_C*(Da+r)**2) + cs*(lcap/(R0*r**2))*np.exp(lcap/(R0*r))/(delta_C*(Da+r));
    n0 = 1.1*np.sqrt(0.5/np.pi) / sigma * np.exp(-1/2 * ((r-0.25 - 1)/sigma)**2);
    n = n0*np.exp(-f_prime*(t-3/(3600*t0)));
    return n


for i in range(t_eval.size):
    print('Time =', t_eval[i])

    h5f = h5py.File('fvm_data/%04d'%i + '.h5', 'r')
    sol = h5f['sol'][:]
    h5f.close()

    h5f = h5py.File('montecarlo_data/%04d'%i + '.h5', 'r')
    R   = h5f['R'][:]
    N   = h5f['N'][:]
    h5f.close()

    #pl.plot(R * R0/1e-9, np.sum(sol[3:-3]) * dr * N, label = 'Particle Model')
    pl.plot(r * R0/1e-9, sol, 'b', label = 'Numerical Simulation');
    #pl.plot(r * R0/1e-9, analytical_1(t_eval[i]), 'k-.', label = 'Travelling Wave Approximation');
    pl.plot(r * R0/1e-9, analytical_2(t_eval[i]), 'r', label = r'Asymptotic Approximation');
    
    pl.title('Time: %3.3f'%(t_eval[i] * t0/3600) + ' hours', fontsize=17)
    pl.xlabel(r'r (in nm)', fontsize=11)
    pl.ylabel(r'N(r,t)', fontsize=11)
    pl.xlim(0.4 * R0/1e-9, 2.5 * R0/1e-9)
    pl.ylim(0, 8)
    pl.legend(fontsize = 11)
    pl.tick_params(axis='both', which='major', labelsize=11);
    pl.savefig('comparison_images/%04d'%i + '.png', bbox_inches = 'tight', dpi=300)
    pl.clf()

