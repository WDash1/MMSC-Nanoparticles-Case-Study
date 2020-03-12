import numpy as NP
from DataVisualiser import DataVisualiser
from PDEModelSimulator import PDEModelSimulator

#The parameter settings that we wish to use for our simulation.
r_min = 0.5;
r_max = 5;

l_cap  = 6e-9;
c_inf_0 = 55.33;

cs   = 5.53e-2;
Vm   = 3.29e-5;
D    = 3.01e-18;
k    = 7.97e-10;
N0   = 8.04e21;
R0   = 3e-9;

end_time = 50;


#This function specifies the starting distribution of our Nanoparticle model.
sigma  = 1/14;
n0 = lambda x: NP.sqrt(0.5/NP.pi) / sigma * NP.exp(-1/2 * ((x - 1)/sigma)**2);


#The number of different injection functions that we wish to use in our
#simualtions. This number should be a value between 0 and 6 (inclusive).
number_of_functions = 6;


#The directory that we wish to write the resulting graph images to.
output_folder = "output_images/";


def injectionFunction1(t):
    # @brief    This example corresponds to the adding of solute to the
    #           solution instantaneously at time t=500.
    # @param t  This parameter corresponds to the current non-dimensionalised
    #           version of time in the simulation that we wish to calculate
    #           current $c_{\infty}(0)$ jump term for.
    # @return   This function must return a number of type double that is >0,
    #           which corresponds to the increase in concentration 
    #           $c_{\infty}(0)$ caused by externally adding more solute to the
    #           system.

    #if(t>5000):
    #    if(t>10000):
    #        if(t>15000):
    #            return 30;
    #        else:
    #            return 20;
    #    else:
    #        return 10;
    #else:
    #    return 0;
    
    #if(t<15000):
    #    return 30*NP.exp(t-15000)
    #else:
    #    return 30;


    if(t<15000):
        return 0*(t-1000)/15000
    else:
        return 0;
def injectionFunction2(t):
    # @brief    This example corresponds to the adding of solute to the
    #           solution instantaneously at time t=5000.
    # @param t  This parameter corresponds to the current non-dimensionalised
    #           version of time in the simulation that we wish to calculate
    #           current $c_{\infty}(0)$ jump term for.
    # @return   This function must return a number of type double that is >0,
    #           which corresponds to the increase in concentration 
    #           $c_{\infty}(0)$ caused by externally adding more solute to the
    #           system.
    #if(t>5000):
    #    if(t>15000):
    #        return 30;
    #    else:
    #        return 15;
    #else:
    #    return 0;
    
    
#    if(t<15000):
#        return 30*NP.log(1+(t*1.71/15000));
#    else:
#        return 30;

    if(t<1000):
        return 0;
    elif(t<16000):
        return 200*(t-1000)/15000
    else:
        return 200;
    

def injectionFunction3(t):
    # @brief    This example corresponds to the adding of solute to the
    #           solution instantaneously at time t=20000.
    # @param t  This parameter corresponds to the current non-dimensionalised
    #           version of time in the simulation that we wish to calculate
    #           current $c_{\infty}(0)$ jump term for.
    # @return   This function must return a number of type double that is >0,
    #           which corresponds to the increase in concentration 
    #           $c_{\infty}(0)$ caused by externally adding more solute to the
    #           system.
    
    #if(t<15000):
    #    return 15*NP.sin(t/1000)+t/1000;
    #else:
    #    return 15*NP.sin(15)+15;
    
    
    if(t<1000):
        return 0;
    elif(t<16000):
        return 400*(t-1000)/15000
    else:
        return 400;

    
    
def injectionFunction4(t):
    # @brief    This example corresponds to the gradual addition of solute 
    #           to the solution using a sigmoid function, starting at time
    #           t=500.
    # @param t  This parameter corresponds to the current non-dimensionalised
    #           version of time in the simulation that we wish to calculate
    #           current $c_{\infty}(0)$ jump term for.
    # @return   This function must return a number of type double that is >0,
    #           which corresponds to the increase in concentration 
    #           $c_{\infty}(0)$ caused by externally adding more solute to the
    #           system.

    #if(t>5000):
#    return 30/(1+NP.exp(-0.00005*(t-10000)))
    #else:
    #    return 0;
    
    
    if(t<1000):
        return 0;
    elif(t<16000):
        return 600*(t-1000)/15000
    else:
        return 600;


def injectionFunction5(t):
    # @brief    This example corresponds to the gradual addition of solute 
    #           to the solution using a sigmoid function, starting at time
    #           t=5000.
    # @param t  This parameter corresponds to the current non-dimensionalised
    #           version of time in the simulation that we wish to calculate
    #           current $c_{\infty}(0)$ jump term for.
    # @return   This function must return a number of type double that is >0,
    #           which corresponds to the increase in concentration 
    #           $c_{\infty}(0)$ caused by externally adding more solute to the
    #           system.
    #if(t>50000):
#    return 15/(1+NP.exp(-0.00005*(t-5000)))+15/(1+NP.exp(-0.00005*(t-15000)));#else:
     #   return 0;
    if(t<1000):
        return 0;
    elif(t<16000):
        return 800*(t-1000)/15000
    else:
        return 800;

def injectionFunction6(t):
    # @brief    This example corresponds to the gradual addition of solute 
    #           to the solution using a sigmoid function, starting at time
    #           t=2000.
    # @param t  This parameter corresponds to the current non-dimensionalised
    #           version of time in the simulation that we wish to calculate
    #           current $c_{\infty}(0)$ jump term for.
    # @return   This function must return a number of type double that is >0,
    #           which corresponds to the increase in concentration 
    #           $c_{\infty}(0)$ caused by externally adding more solute to the
    #           system.if(t>20000):
    #if(t>200000):
#    return 30/(1+NP.exp(-0.00005*(t-10000)));
    if(t<1000):
        return 0;
    elif(t<33000):
        return 1000*(t-1000)/15000
    else:
        return 1000;

    #else:
    #    return 0;


def linearInjectionFunction(t, start_time, stop_time, min_value, max_value):
    if(t<=start_time):
        return min_value;
    elif(t<stop_time):
        return min_value+((t-start_time)/(stop_time-start_time)*(max_value-min_value));
    else:
        return max_value;
    

#This array specifies the description strings in the legend for each line on
#the output graphs.
line_key_strings = ["0",
                    "200",
                    "400",
                    "600",
                    "800",
                    "1000"];

#This array specifies the colour of each line on the output graphs.
line_colours = ["red",
                "green",
                "blue",
                "yellow",
                "purple",
                "orange"];

                
injection_start_time=2000;
injection_stop_time=64000;

injection_amounts=[0,200,400,600,800,1000];

                
injection_functions = [lambda t: linearInjectionFunction(t, injection_start_time, injection_stop_time,0,injection_amounts[0]),
                       lambda t: linearInjectionFunction(t, injection_start_time, injection_stop_time,0,injection_amounts[1]),
                lambda t: linearInjectionFunction(t, injection_start_time, injection_stop_time,0,injection_amounts[2]),
                lambda t: linearInjectionFunction(t, injection_start_time, injection_stop_time,0,injection_amounts[3]),
                lambda t: linearInjectionFunction(t, injection_start_time, injection_stop_time,0,injection_amounts[4]),
                lambda t: linearInjectionFunction(t, injection_start_time, injection_stop_time,0,injection_amounts[5])];
                
#This array stores referances to each of the injection functions we wish to
#use in simulations of our PDE model.                    
#injection_functions = [injectionFunction1, 
#                       injectionFunction2, 
#                       injectionFunction3, 
#                       injectionFunction4, 
#                       injectionFunction5, 
#                       injectionFunction6];
                       
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
visualiser = DataVisualiser(1,   0.4 * R0/1e-9, 4.0 * R0/1e-9, 'Nanoparticle Radius (in nanometres)',
                                 0, 12, 'N(r,t)', 
                                 line_key_strings[0:number_of_functions],
                                 line_colours[0:number_of_functions]);
                          
                            
                             
#For each time step create a graph containing the N distribution produced using
#each of the different injection functions and export them to an external
#image file.

mean_values = NP.empty(number_of_functions, dtype=object);
variance_values = NP.empty(number_of_functions, dtype=object);                 

moment_calculator = lambda n_data, moment_number: NP.dot(
                n_data, NP.power(r_values, moment_number))/sum(n_data);


mean_calculator = lambda n_data: moment_calculator(n_data, 1);

variance_calculator = lambda n_data: moment_calculator(n_data, 2) - mean_calculator(n_data)**2;





max_means=list(range(0,number_of_functions));
min_means=list(range(0,number_of_functions));

max_vars=list(range(0,number_of_functions));
min_vars=list(range(0,number_of_functions));
                                 

for i in range(number_of_functions):
    n_values = list(map(lambda j: solution_data[i][:,j], 
                        list(range(0,time_values.size))));
    
                        
    mean_values[i] = list(map(lambda j: mean_calculator(solution_data[i][:,j]), 
                        list(range(0,time_values.size))));
               
    max_means[i] = max(mean_values[i]);
    min_means[i] = min(mean_values[i]);           
               
    variance_values[i] = list(map(lambda j: variance_calculator(solution_data[i][:,j]), 
                        list(range(0,time_values.size))));
    max_vars[i] = max(variance_values[i]);
    min_vars[i] = min(variance_values[i]);
    
mean_visualiser = DataVisualiser(1,   min(time_values), max(time_values), 'Time Since Nucleation (in seconds)',
                                 min(min_means)*0.9, max(max_means)*1.1, 'Average Nanoparticle Radius (in nanometres)', 
                                 line_key_strings[0:number_of_functions],
                                 line_colours[0:number_of_functions]);
                                 
                                 
variance_visualiser = DataVisualiser(1,   min(time_values), max(time_values), 'Time Since Nucleation (in seconds)',
                                 min(min_vars)*0.9, max(max_vars)*1.1, 'Variance of Nanoparticle Radius (in nanometres)', 
                                 line_key_strings[0:number_of_functions],
                                 line_colours[0:number_of_functions]);
                   
#add_dataset_function = lambda dataset: map(mean_visualiser.addData(time_values, dataset[i],i), list(range(0,dataset.size)));


#add_dataset_function(mean_values);

for i in range(number_of_functions):
    mean_visualiser.addData(time_values, mean_values[i],i);
    variance_visualiser.addData(time_values, variance_values[i],i)

mean_visualiser.exportGraph(' ', 'mean_graph.png', True);
variance_visualiser.exportGraph(' ', 'variance_graph.png', True);
mean_visualiser.clearData();
variance_visualiser.clearData();


               
#    mean_values[i] = list(map(lambda j: mean_calculator(n_values[j,:]), list(range(0,time_values.size))));
 #   variance_values[i] = list(map(lambda j: variance_calculator(n_values[j,:]), list(range(0,time_values.size))));
          
for i in range(time_values.size):
    print('Time =', time_values[i])

    for j in range(number_of_functions):
        current_N_data = solution_data[j][:,i];    
        visualiser.addData(1/1e-9*r_values, current_N_data, j);   

    visualiser.exportGraph('Time Since Nucleation: %3.3f'%(time_values[i]) + ' seconds', output_folder+'/%04d'%i + '.png', False);
    visualiser.clearData();
    