##  @brief  This file defines the PDEModelSimulator class.


import numpy as NP;
from scipy.integrate import solve_ivp

from Simulation_Utilities import riemannSolver;
from Simulation_Utilities import reconstructMinmod;


##  @brief  This class provides a way of initialising, running and retriving
#           the corresponding output of numerical simulations of the PDE model
#           for our system of nanoparticles.
class PDEModelSimulator:
    
    ##  @brief          The constructor for this class. This method is 
    #                   responsibleb for initialising all the necessary
    #                   properties of this class in order for simultations 
    #                   of our model to be executed. However it is important
    #                   to note that this method itself does compute any such
    #                   simulations.
    #   @param t_final  This parameter must be a number of type double which
    #                   is >0 and corresponds to the final
    #                   (non-dimensionalised) time at which
    #                   numerical simulations of our PDE model are terminated.
    #   @param r_min    This parameter must be a number of type double which
    #                   corresponds to the theshold size at which nanoparticles
    #                   are considered to have discolved back into the
    #                   bulk solution. It is important to note that this
    #                   parameter is relative to the initial nanoparticle size
    #                   (denoted by the parameter "R0" ). In addition, it is
    #                   important to note that this parameter must be strictly
    #                   less than the parameter "r_max".
    #   @param r_max    This parameter must be a number of type double, which
    #                   corresponds to the maximum size of particles that this
    #                   simulation will consider. It is important to note
    #                   that this parameter must be strictly greater than 
    #                   r_min and must be sufficiently large so that no
    #                   nanoparticles in the simulation every actually attain
    #                   this value. In addition, it is important to note
    #                   that this parameter is relative to the parameter "R0".
    #   @param l_cap    This parameter must be a positive number of type double
    #                   which corresponds to the capillarly length coefficient
    #                   in our model. (For referance, this parameter should 
    #                   typically be of the order 10^-5).
    #   @param c_inf_0  This parameter must be a positive number of type double
    #                   that dictates the bulk concentraiton of the solute
    #                   in the solution at time t=0.
    #   @param cs       This parameter must be a positive number of type double
    #                   which corresponds to the solubility of solute and is
    #                   used in the Gibbs-Thompson equation.
    #   @param Vm       This parameter must be a positive number of type double
    #                   which dictates the molar volume (that is the volume
    #                   occupied by 1 mol of the solute species).
    #   @param D        This parameter must be a positive number of type double
    #                   which corresponds to the diffusion coefficient
    #                   associated with our solute species, when it is
    #                   diffusion through our solution.
    #   @param k        This parameter must be a positive number of type double
    #                   which corresponds to the rate of reactivity of the
    #                   chemcial reaction which is taking place on the surface
    #                   of each of the nanoparticles.
    #   @param N0       This parameter must be a positive number of type double
    #                   which corresponds to the total volume of the bath in
    #                   which the nanoparticle reaction process is taking
    #                   place.
    #   @param R0       This parameter must be a positive integer of type
    #                   double that corresponds to the average initial size
    #                   of the nanoparticles in our distribution at time t=0.
    #                   This parameter is used in order to non-dimensionalise
    #                   the lengthscale of our partial differential equation
    #                   and should typically be of the order of 10^-9.
    def __init__(self, t_final, r_min, r_max, 
                             l_cap, c_inf_0, cs, Vm, D, k, N0, R0):

        #Interally store all of the provided parameters and calculate any
        #additional parameters that are required for simulations.
        self.__r_min = r_min;
        self.__r_max = r_max;
        self.__N_r   = 1024;
        self.__N_g   = 3;
        self.__dr    = (self.__r_max - self.__r_min) / self.__N_r;
        self.__r     = self.__r_min + (0.5 + NP.arange(-self.__N_g, self.__N_r 
                                                   + self.__N_g)) * self.__dr;

        self.__lcap  = l_cap;
        self.__cinf0 = c_inf_0;
        self.__cinfu = self.__cinf0;
        self.__cs   = cs;
        self.__Vm   = Vm;
        self.__D    = D;
        self.__k    = k;
        self.__N0   = N0;
        self.__R0   = R0;
        self.__beta = 4 * NP.pi * self.__N0 / (3 * self.__Vm);

        self.__Da      = self.__D/(self.__k * self.__R0 );
        self.__delta_C = self.__cinf0 - self.__cs * NP.exp(self.__lcap/(self.__R0));
        self.__t0      = (self.__R0)**2 / (self.__Vm * self.__D * self.__delta_C);
        self.__t_final = t_final;
        self.__t_eval  = NP.append(NP.linspace(0, 5, 51), 
                                 NP.arange(6, self.__t_final + 1));

    
    ##  @brief                      This method is responsible for executing a
    #                               simulation of the PDE model using the
    #                               settings given to the constructor and 
    #                               returning the numerical result.
    #   @param n0                   This parameter should be a function handle
    #                               that corresponds to the initial
    #                               distribution of particle sizes at time t=0.
    #   @param injection_function    This parameter should be a function handle
    #                               that takes a dimensionalised time argument
    #                               and returns a number of type double, which
    #                               is strictly positive. This parameter
    #                               corresponds to the amount of additional
    #                               solute that has been externally added to
    #                               the solution throughout the course of the
    #                               simulation of the model.
    #   @return                     This method returns a two dimensional
    #                               array of size AxB, where A corresponds to
    #                               the amount of different r values used at
    #                               each time step of the simulations of this
    #                               method (and can be accessed via the
    #                               getRValues method of this class) and B
    #                               corresponds to the amount of different t
    #                               values used for simulations of this method
    #                               (and can be accessed via the getTValues
    #                               method of this class).
    def simulateModel(self, n0, injection_function):
        n0_data_values = list(map(n0, self.__r));
        self.__cinfu = self.__cinf0;
        df_dt_function = lambda t,f: self.__calculateDfDt(t, f, injection_function);
        simulation_output = solve_ivp(df_dt_function, (0, self.__t_final), n0_data_values, t_eval=self.__t_eval);
    
        return simulation_output.y;


    ##  @brief          This method can be used to retrieve the dimensionalised
    #                   r values that our simulation mesh grid is comprised of.
    #   @return         This method returns a one dimensional array of type
    #                   double that is comprised on the ordered radius values
    #                   used in the mesh of our domain for which numerical
    #                   solutions have been produced by this class.
    def getRValues(self):
        r_values = self.__R0*self.__r;
        return r_values;

    ##  @brief          This method can be used to retrieved the
    #                   the dimensionalised t values that our simulation mesh
    #                   grid is comprised of.
    #   @return         This method returns a one dimensional array of type
    #                   double that is comprosed of the time values used
    #                   at each step of our numerical approximation of 
    #                   the solution of our PDE model.
    def getTValues(self):
        t_values = self.__t0*self.__t_eval;
        return t_values;
    
    ##  @brief                      This method is used internally by the class
    #                               in order to calculate the derivative with
    #                               respect to time of our distribution N for
    #                               a given r and t value.
    #   @param f                    This parameter is a one dimensional array
    #                               of type double that corresponds to the
    #                               values of the distribution N at time step t
    #                               in the simulation.
    #   @param t                    This parameter corresponds to the current
    #                               non-dimensionalised time value that we are
    #                               evaluating the partial derivative of N
    #                               with respect to t at in our numerical
    #                               simulation of our PDE model.
    #   @param injection_function   This parameter must be a function handle
    #                               that takes the current dimensionalised time
    #                               value as an input and returns the 
    #                               total amount of solute that have been added
    #                               to to the bulk solution before the given
    #                               time.
    def __calculateDfDt(self, t, f, injection_function):
        self.__cinfu  += self.__beta * self.__R0**3 * NP.sum(f[:self.__N_g] * self.__r[:self.__N_g]**3) * self.__dr;
        f[:self.__N_g] = 0;
        if(NP.sum(f[-self.__N_g:] * self.__r[-self.__N_g:]**3) * self.__dr > 1e-8):
            raise Exception('Increase r_max!')

        f_left_plus_eps, f_right_minus_eps = reconstructMinmod(f)
        f_left_minus_eps = NP.roll(f_right_minus_eps, 1)

        cinf     = self.__cinfu - self.__beta * self.__R0**3 * NP.sum(f[self.__N_g:-self.__N_g] * self.__r[self.__N_g:-self.__N_g]**3) * self.__dr
    
        #Modify the boundary condition in accordance with any additional solute
        #has been externally added to the system.
        cinf = cinf + injection_function(t*self.__t0);

        velocity = (cinf - self.__cs * NP.exp(self.__lcap / (self.__R0 * self.__r))) / (self.__delta_C * (self.__Da + self.__r))
        f_left   = riemannSolver(f_left_minus_eps, 
                              f_left_plus_eps, 
                              velocity
                             )
    
        left_flux  = velocity * f_left
        right_flux = NP.roll(left_flux, -1)
    
        df_dt = -(right_flux - left_flux) / self.__dr
        return df_dt

    

