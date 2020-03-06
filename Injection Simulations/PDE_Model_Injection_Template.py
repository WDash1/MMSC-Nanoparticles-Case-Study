import numpy as NP
from DataVisualiser import DataVisualiser
from PDEModelSimulator import PDEModelSimulator

#The parameter settings that we wish to use for our simulation.
r_min = 0.5;
r_max = 3;

l_cap  = 6e-9;
c_inf_0 = 55.33;

cs   = 5.53e-2;
Vm   = 3.29e-5;
D    = 3.01e-18;
k    = 7.97e-10;
N0   = 8.04e21;
R0   = 3e-9;

end_time = 200;


#This function specifies the starting distribution of our Nanoparticle model.
sigma  = 1/14;
n0 = lambda x: NP.sqrt(0.5/NP.pi) / sigma * NP.exp(-1/2 * ((x - 1)/sigma)**2);


#The number of different injection functions that we wish to use in our
#simualtions. This number should be a value between 0 and 6 (inclusive).
number_of_functions = 2;


#The directory that we wish to write the resulting graph images to.
output_folder = "output_images/";


def injectionFunction1(t):
    # @TODO     Modify this function to find the optimal way to add more
    #           solute to the solution.
    # @brief    This example simply uses the Heaviside function which
    #           corresponds to an increase in concentration of 5 mols per metre
    #           of solute that occurs instantaneously at dimensionalised
    #           time t=10.
    # @param t  This parameter corresponds to the current non-dimensionalised
    #           version of time in the simulation that we wish to calculate
    #           current $c_{\infty}(0)$ jump term for.
    # @return   This function must return a number of type double that is >0,
    #           which corresponds to the increase in concentration 
    #           $c_{\infty}(0)$ caused by externally adding more solute to the
    #           system.

    if(t>10):
        return 5;
    else:
        return 0;

def injectionFunction2(t):
    # @TODO     Modify this function to find the optimal way to add more 
    #           solute to the solution.
    #           This example adds more solute to the solution, so that the
    #           far field concentration increases linearly with time from t=5 
    #           to t=10.
    # @param t  This parameter corresponds to the current non-dimensionalised
    #           version of time in the simulation that we wish to calculate
    #           current $c_{\infty}(0)$ jump term for.
    # @return   This function must return a number of type double that is >0,
    #           which corresponds to the increase in concentration 
    #           $c_{\infty}(0)$ caused by externally adding more solute to the
    #           system.
    if(t>5):
        if(t<10):
            return (t-5)*2;
        else:
            return (10-5)*2;
    else:
        return 0;

def injectionFunction3(t):
    # @TODO     Modify this function to find the optimal way to add more 
    #           solute to the solution.
    return 0;

def injectionFunction4(t):
    # @TODO     Modify this function to find the optimal way to add more 
    #           solute to the solution.
    return 0;

def injectionFunction5(t):
    # @TODO     Modify this function to find the optimal way to add more 
    #           solute to the solution.
    return 0;

def injectionFunction6(t):
    # @TODO     Modify this function to find the optimal way to add more 
    #           solute to the solution.
    return 0;

#This array specifies the description strings in the legend for each line on
#the output graphs.
line_key_strings = ["Injection function 1 description",
                    "Injection function 2 description",
                    "Injection function 3 description",
                    "Injection function 4 description",
                    "Injection function 5 description",
                    "Injection function 6 description"];

#This array specifies the colour of each line on the output graphs.
line_colours = ["red",
                "green",
                "blue",
                "yellow",
                "purple",
                "orange"];

#This array stores referances to each of the injection functions we wish to
#use in simulations of our PDE model.                    
injection_functions = [injectionFunction1, 
                       injectionFunction2, 
                       injectionFunction3, 
                       injectionFunction4, 
                       injectionFunction5, 
                       injectionFunction6];
                       
#Setup a simulation environement and retrieve the output r and t values that
#that will be used by any such simulation.
simulation_environment = PDEModelSimulator( end_time, r_min, r_max, 
                                    l_cap, c_inf_0, cs, Vm, D, k, N0, R0)
                        
time_values = simulation_environment.getTValues();
r_values = simulation_environment.getRValues();

                       
#Create a function that performs numerical simulations of the PDE model for
#a given injection function index.                       
numerical_simulator = lambda index: simulation_environment.simulateModel(n0, 
                                    injection_functions[index]);

#Perform simulations of the PDE model for each different injection function.
solution_data = list(map(numerical_simulator, 
                           NP.arange(number_of_functions)));

#Create a DataVisualiser object instance in order to export the numeric
#results into graphical form.
visualiser = DataVisualiser(1,   0.4 * R0/1e-9, 2.5 * R0/1e-9, 'r(in nm)',
                                 0, 6, 'N(r,t)', 
                                 line_key_strings[0:number_of_functions],
                                 line_colours[0:number_of_functions]);

                            
                             
#For each time step create a graph containing the N distribution produced using
#each of the different injection functions and export them to an external
#image file.
for i in range(time_values.size):
    print('Time =', time_values[i])

    for j in range(number_of_functions):
        current_N_data = solution_data[j][:,i];    
        visualiser.addData(1/1e-9*r_values, current_N_data, j);   

    visualiser.exportGraph('Time = %3.3f'%(time_values[i]) + ' seconds', output_folder+'/%04d'%i + '.png');
    visualiser.clearData();
    