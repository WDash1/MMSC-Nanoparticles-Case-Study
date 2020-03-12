from params import *
import h5py
from scipy.special import erfinv
from scipy.integrate import solve_ivp
import pylab as pl

R   = r_mean + sigma * np.sqrt(2) * erfinv(2 * np.random.rand(num) - 1)
def dR_dt(t, R):
    # Identifies particles which are less than 0.5
    ind     = R < r_min
    nind    = np.invert(ind)
    # N       = np.where(nind)[0].size
    # if(num != N):
    #    print(cinf0 + beta * np.sum(((R0 * R)[ind])**3) / (num - N))
    csolute = cinf0 - beta * np.sum(((R0 * R)[nind])**3) / num
    dRdt    = (csolute - cs * np.exp(lcap / (R0 * R)))/(delta_C * (Da + R))
    
    dRdt[ind] = 0
    return dRdt

sol = solve_ivp(dR_dt, (0, t_final), R, t_eval=t_eval)
for i in range(sol.y.shape[1]):
    solution = sol.y[:, i]
    solution = solution[solution > r_min]
    N, R_bin = np.histogram(solution, 100, density = True)
    h5f = h5py.File('montecarlo_data/%04d'%i + '.h5', 'w')
    h5f.create_dataset('R', data = (R_bin[1:] + R_bin[:-1])/2)
    h5f.create_dataset('N', data = N)
    h5f.close()
