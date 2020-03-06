import numpy as NP;
from scipy.integrate import solve_ivp

from Simulation_Utilities import riemannSolver;
from Simulation_Utilities import reconstructMinmod;


class PDEModelSimulator:
    
    def __init__(self, t_final, r_min, r_max, 
                             l_cap, c_inf_0, cs, Vm, D, k, N0, R0):

        self.r_min = r_min;
        self.r_max = r_max;
        self.N_r   = 1024;
        self.N_g   = 3;
        self.dr    = (self.r_max - self.r_min) / self.N_r;
        self.r     = self.r_min + (0.5 + NP.arange(-self.N_g, self.N_r + self.N_g)) * self.dr;

        self.lcap  = l_cap;
        self.cinf0 = c_inf_0;
        self.cinfu = self.cinf0;
        self.cs   = cs;
        self.Vm   = Vm;
        self.D    = D;
        self.k    = k;
        self.N0   = N0;
        self.R0   = R0;
        self.beta = 4 * NP.pi * self.N0 / (3 * self.Vm);

        self.Da      = self.D/(self.k * self.R0 );
        self.delta_C = self.cinf0 - self.cs * NP.exp(self.lcap/(self.R0));
        self.t0      = (self.R0)**2 / (self.Vm * self.D * self.delta_C);
        self.t_final = t_final;
        self.t_eval  = NP.append(NP.linspace(0, 5, 51), NP.arange(6, self.t_final + 1));

    def simulateModel(self, n0, injection_function):
        n0_data_values = list(map(n0, self.r));
        self.cinfu = self.cinf0;
        df_dt_function = lambda t,f: self.calculateDfDt(t, f, injection_function);
        simulation_output = solve_ivp(df_dt_function, (0, self.t_final), n0_data_values, t_eval=self.t_eval);
    
        return simulation_output.y;


    def getRValues(self):
        r_values = self.R0*self.r;
        return r_values;

    def getTValues(self):
        t_values = self.t0*self.t_eval;
        return t_values;
    
    def calculateDfDt(self, t, f, injection_function):
        self.cinfu  += self.beta * self.R0**3 * NP.sum(f[:self.N_g] * self.r[:self.N_g]**3) * self.dr;
        f[:self.N_g] = 0;
        if(NP.sum(f[-self.N_g:] * self.r[-self.N_g:]**3) * self.dr > 1e-8):
            raise Exception('Increase r_max!')

        f_left_plus_eps, f_right_minus_eps = reconstructMinmod(f)
        f_left_minus_eps = NP.roll(f_right_minus_eps, 1)

        cinf     = self.cinfu - self.beta * self.R0**3 * NP.sum(f[self.N_g:-self.N_g] * self.r[self.N_g:-self.N_g]**3) * self.dr
    
        #Modify the boundary condition in accordance with any additional solute
        #has been externally added to the system.
        cinf = cinf + injection_function(t*self.t0);

        velocity = (cinf - self.cs * NP.exp(self.lcap / (self.R0 * self.r))) / (self.delta_C * (self.Da + self.r))
        f_left   = riemannSolver(f_left_minus_eps, 
                              f_left_plus_eps, 
                              velocity
                             )
    
        left_flux  = velocity * f_left
        right_flux = NP.roll(left_flux, -1)
    
        df_dt = -(right_flux - left_flux) / self.dr
        return df_dt

    

