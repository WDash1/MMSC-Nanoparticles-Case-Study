import h5py
from scipy.special import erfinv
from scipy.integrate import solve_ivp
import pylab as pl

import numpy as NP;

class RModelSimulator:
    def __init__(self, t_final, r_mean, l_cap, c_inf_0, cs, Vm, D, k, N0, R0):
        self.t_final = t_final;
        self.r_mean = r_mean;
        self.l_cap = l_cap;
        self.c_inf_0 = c_inf_0;
        self.cs = cs;
        self.Vm = Vm;
        self.D = D;
        self.k = k;
        self.N0 = N0;
        self.R0 = R0;
        
        self.Da      = self.D/(self.k * self.R0 );
        self.delta_C = self.cinf0 - self.cs * NP.exp(self.lcap/(self.R0));
        self.t0      = (self.R0)**2 / (self.Vm * self.D * self.delta_C);
        self.t_eval  = NP.append(NP.linspace(0, 5, 51), 
                                 NP.arange(6, self.t_final + 1));
        self.t_final = t_final;
        


    def calculateDRDt(self, t, R, injection_function):
        # Identifies particles which are less than 0.5
        ind     = R < self.r_min;
        nind    = NP.invert(ind);
        # N       = np.where(nind)[0].size
        # if(num != N):
        #    print(cinf0 + beta * np.sum(((R0 * R)[ind])**3) / (num - N))
        csolute = self.cinf0 - self.beta * NP.sum(((self.R0 * R)[nind])**3) / self.num;
        
        #Apply any excess solute added to the solution.
        csolute = csolute + injection_function(t*self.t0);
        
        dRdt    = (csolute - self.cs * NP.exp(self.lcap / (self.R0 * self.R)))/(self.delta_C * (self.Da + R));
    
        dRdt[ind] = 0;
        return dRdt;
    
    def simulateModel(self, n0, injection_function):

        R   = r_mean + sigma * np.sqrt(2) * erfinv(2 * np.random.rand(num) - 1)
        
        n0_data_values = list(map(n0, self.r));

        dr_dt_function = lambda t,r: self.calculateDRDt(t, r, injection_function);
        
        sol = solve_ivp(dr_dt_function, (0, t_final), R, t_eval=t_eval);

        
        simulation_output = solve_ivp(df_dt_function, (0, self.t_final), n0_data_values, t_eval=self.t_eval);
    
        return simulation_output.y;

sol = solve_ivp(dR_dt, (0, t_final), R, t_eval=t_eval)
for i in range(sol.y.shape[1]):
    solution = sol.y[:, i]
    solution = solution[solution > r_min]
    N, R_bin = np.histogram(solution, 100, density = True)
    h5f = h5py.File('montecarlo_data/%04d'%i + '.h5', 'w')
    h5f.create_dataset('R', data = (R_bin[1:] + R_bin[:-1])/2)
    h5f.create_dataset('N', data = N)
    h5f.close()
