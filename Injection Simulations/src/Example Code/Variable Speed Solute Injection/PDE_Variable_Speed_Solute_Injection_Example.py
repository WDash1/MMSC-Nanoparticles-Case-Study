##  @brief      This file includes the source code that was used to generate
#               the distribution images for solute injections of different
#               speeds, whose graphs are shown in slide 33 of the
#               presentation.

import numpy as NP

import sys

sys.path.append('../../Utilities/')

from PDEModelSimulator import PDEModelSimulator
from DistributionDataProcessor import DistributionDataProcessor;
from DistributionDataVisualiser import DistributionDataVisualiser;

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

end_time = 7;


#This function specifies the starting distribution of our Nanoparticle model.
sigma  = 1/14;
n0 = lambda x: NP.sqrt(0.5/NP.pi) / sigma * NP.exp(-1/2 * ((x - 1)/sigma)**2);


#The number of different injection functions that we wish to use in our
#simualtions. This number should be a value between 0 and 6 (inclusive).
number_of_functions = 4;


#The directory that we wish to write the resulting graph images to.
output_folder = "Distribution_Plots/";

#The size of the fonts to be used for the axis and legend labels of the output
#graphs.
label_font_size = 13;

#The size of the fonts to be used for the title of the output graphs.
title_font_size = 17;


#This function corresponds to no solute injection.
def injectionFunction1(t):
    return 0;

#This function corresponds to a linear solute inejection.
def injectionFunction2(t, start_time, end_time, injection_amount):
    if(t>start_time):
        if(t>=end_time):
            return injection_amount;
        else:
            return (t-start_time)/(end_time-start_time)*injection_amount;
    else:
        return 0;
    
#This function corresponds to a quadratic solute injection.
def injectionFunction3(t, start_time, end_time, injection_amount):
    if(t>start_time):
        if(t>=end_time):
            return injection_amount;
        else:
            return (((t-start_time)/(end_time-start_time))**2)*injection_amount;
    else:
        return 0;
    
#This function corresponds to a square root solute injection.
def injectionFunction4(t, start_time, end_time, injection_amount):
    if(t>start_time):
        if(t>=end_time):
            return injection_amount;
        else:
            return NP.sqrt(t-start_time)/(NP.sqrt(end_time-start_time))*injection_amount;
    else:
        return 0;
        
    
#This function corresponds to a logarithmic solute injection.
def injectionFunction5(t, start_time, end_time, injection_amount):
    if(t>start_time):
        if(t>=end_time):
            return injection_amount;
        else:
            return NP.log(t+1-start_time)/NP.log(end_time+1-start_time)*(injection_amount);
    else:
        return 0;


#This function corresponds to an exponential solute injection.    
def injectionFunction6(t, start_time, end_time, injection_amount):
    if(t>start_time):
        if(t>=end_time):
            return injection_amount;
        else:
            scale = NP.log(0.001);
            return NP.exp(scale*(end_time-t)/(end_time-start_time))*injection_amount;
    else:
        return 0;


    

#This array specifies the description strings in the legend for each line on
#the output graphs.
line_key_strings = ["original",
                    "linear",
                    "quadratic",
                    "square root",
                    "exponential",
                    "logarithmic"];

#This array specifies the colour of each line on the output graphs.
line_colours = ["black",
                "red",
                "green",
                "blue",
                "yellow",
                "purple"];

#This array specifies the type of each line on the output graphs.
line_styles = ["--",
               "--",
               "--",
               "--",
               "--",
               "--"];

              
#This parameter dictates the dimensionalised time (in seconds) that we will
#start and stop injecting solute (respectively).
injection_start_time=2000;
injection_stop_time=5000;

#This parameter dictates the change in concentration $c_{\infty(0)}$ caused by
#the injection of solute into the reaction.
injection_amount=50;
                

#This array stores referances to each of the injection functions we wish to
#use in simulations of our PDE model.                    
injection_functions = [injectionFunction1, 
                       lambda t: injectionFunction2(t, injection_start_time, injection_stop_time, injection_amount), 
                       lambda t: injectionFunction3(t, injection_start_time, injection_stop_time, injection_amount), 
                       lambda t: injectionFunction4(t, injection_start_time, injection_stop_time, injection_amount), 
                       lambda t: injectionFunction5(t, injection_start_time, injection_stop_time, injection_amount), 
                       lambda t: injectionFunction6(t, injection_start_time, injection_stop_time, injection_amount)];
                       
#Setup a simulation environement and retrieve the output r and t values that
#that will be used by any such simulation.
simulation_environment = PDEModelSimulator( end_time, r_min, r_max, 
                                    l_cap, c_inf_0, cs, Vm, D, k, N0, R0)

#Retrieve the time and r values used in the simulation and scale them
#appropriately.
time_values = simulation_environment.getTValues()/3600;
r_values = simulation_environment.getRValues()*1e9;

          
#Put the data into a processor so that we can create graphs an calculate
#some useful statistics.
data_processor = DistributionDataProcessor(r_values, time_values, 
                                           number_of_functions); 
for j in range(number_of_functions):
    data = simulation_environment.simulateModel(n0, injection_functions[j]);
    data_processor.addData(data, j);

#Create a visualiser for the simulation data.
data_visualiser = DistributionDataVisualiser("r (in nm)", "N(r,t)", 
                                             line_key_strings[0:number_of_functions],
                                             line_colours[0:number_of_functions],
                                             line_styles[0:number_of_functions],
                                             label_font_size,
                                             title_font_size,
                                             data_processor);   

                                             
#Generate a graph of how the mean of the distribution evolves with time.
data_visualiser.exportMeanTimeGraph(1, "", "t (in hours)", "Average r (in nm)", 
                                    "./Statistics/mean_graph.png", True,
                                    True, 0, 0);

#Generate a graph of how the variance of the distirbution evolves with time.
data_visualiser.exportVarianceTimeGraph(2, "", "t (in hours)", 
                                        "Variance of r (in nm$^2$)",
                                        "./Statistics/variance_graph.png",
                                        True, True, 0, 0);

#Generate graphs of the distribution at each time step in the simulation.                                        
data_visualiser.exportDistributionGraphAll(output_folder, False, False, 0, 7,
                                           True, 0, 0);	