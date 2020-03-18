##  @brief      This file includes the source code that was used to generate
#               the distribution images for solute injection at different times
#               in the reaction, whose graphs are shown in slides 31 and 32 of
#               the presentation.


import numpy as NP

import sys
sys.path.append('../../Utilities/')

from PDEModelSimulator import PDEModelSimulator
from DistributionDataProcessor import DistributionDataProcessor;
from DistributionDataVisualiser import DistributionDataVisualiser;


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

end_time = 600


#The size of the fonts to be used for the axis and legend labels of the output
#graphs.
label_font_size = 13;

#The size of the fonts to be used for the title of the output graphs.
title_font_size = 17;

#This function specifies the starting distribution of our Nanoparticle model.
sigma  = 1/14;
n0 = lambda x: NP.sqrt(0.5/NP.pi) / sigma * NP.exp(-1/2 * ((x - 1)/sigma)**2);


#The number of different injection functions that we wish to use in our
#simualtions. This number should be a value between 0 and 6 (inclusive).
number_of_functions = 6;


#The directory that we wish to write the resulting graph images to.
output_folder = "Raw_Images/";


#A general step function that we will use to rapidly sharply inject solute
#into the system.
def sharpInjectionFunction(t, threshold, start_value, end_value):
    if(t<=threshold):
        return start_value;
    else:
        return end_value;

    

#This array specifies the description strings in the legend for each line on
#the output graphs.
line_key_strings = ["Original", "0.278", "13.889", "27.778", "55.556", "138.889"];
                    

#This array specifies the colour of each line on the output graphs.
line_colours = ["black",
                "red",
                "green",
                "yellow",
                "blue",
                "purple",
                "orange"];

                
#This array specifies the type of each line on the output graphs.
line_styles = ["--",
               "--",
               "--",
               "--",
               "--",
               "--"];

#The total change in c_{\infty}(0) caused by the injection.
injection_amount = 50;

#The time (in dimensionalised seconds at which we inject solute).
injection_times = [1000, 50000, 100000, 200000, 500000];


#Generate all the necessarly solute injection functions using the settings
#above.                
injection_functions = [lambda t: 0,
                       lambda t: sharpInjectionFunction(t, injection_times[0], 0, injection_amount),
                       lambda t: sharpInjectionFunction(t, injection_times[1], 0, injection_amount),
                       lambda t: sharpInjectionFunction(t, injection_times[2], 0, injection_amount),
                       lambda t: sharpInjectionFunction(t, injection_times[3], 0, injection_amount),
                       lambda t: sharpInjectionFunction(t, injection_times[4], 0, injection_amount)];   

                       
#Setup a simulation environement and retrieve the output r and t values that
#that will be used by any such simulation.
simulation_environment = PDEModelSimulator( end_time, r_min, r_max, 
                                    l_cap, c_inf_0, cs, Vm, D, k, N0, R0)
                        
#Retrieve the simulation r and t values and scale appropriately.
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
                                    "./Stats/mean_graph.png", True,
                                    True, 0, 0);

#Generate a graph of how the variance of the distirbution evolves with time.
data_visualiser.exportVarianceTimeGraph(2, "", "t (in hours)", 
                                        "Variance of r (in nm$^2$)",
                                        "./Stats/variance_graph.png", True,
                                        True, 0, 0);

#Generate graphs of the distribution at each time step in the simulation.                                        
data_visualiser.exportDistributionGraphAll(output_folder, False, True, 0, 0);	 