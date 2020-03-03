import numpy as np
import h5py
from scipy.special import erfinv
from scipy.integrate import odeint
import pylab as pl

lcap  = 6e-9
cinf0 = 55.33

cs   = 5.53e-2
Vm   = 3.29e-5
D    = 3.01e-18
k    = 7.97e-10
N0   = 8.04e21
R0   = 3e-9
beta = 4 * np.pi * N0 / (3 * Vm)
Da   = D/(k * R0)

delta_C = cinf0 - cs * np.exp(lcap/R0)
t0      = (R0**2 / (Vm * D * delta_C))

# Number of particles that we are sampling:
num = 10000
R   = 1 + (np.sqrt(2) / 10) * erfinv(2 * np.random.rand(num) - 1)

def dR_dt(R, t):
    dRdt = (cinf0 - beta * np.sum((R0 * R)**3)/R.size - cs * np.exp(lcap / (R0 * R)))/(delta_C * (Da + R))
    return dRdt

N, R_bin = np.histogram(R, 50, density = True)

h5f = h5py.File('montecarlo_data/0000.h5', 'w')
h5f.create_dataset('R', data = (R_bin[1:] + R_bin[:-1])/2)
h5f.create_dataset('N', data = N)
h5f.close()

t_final = 200
for T in np.arange(1, t_final):
    print('Computing For Time =', T)
    N = 100
    if(T>1):
        while(sol[1]['message'] != 'Integration successful.'):
            print('N is now =', N)
            sol = odeint(dR_dt, old_sol, np.linspace(T-1, T, N), full_output = 1)
            N  *= 10

        old_sol = sol[0][-1]
        # Reseting the message:
        sol[1]['message'] = '1'
    else:
        sol = odeint(dR_dt, R, np.linspace(T-1, T, N), full_output = 1)
        assert(sol[1]['message'] == 'Integration successful.')
        old_sol = sol[0][-1]
        # Reseting the message:
        sol[1]['message'] = '1'
    
    N, R_bin = np.histogram(old_sol, 50, density = True)

    h5f = h5py.File('montecarlo_data/%04d'%T + '.h5', 'w')
    h5f.create_dataset('R', data = (R_bin[1:] + R_bin[:-1])/2)
    h5f.create_dataset('N', data = N)
    h5f.close()
