#   @brief  This file acts as an example of how you might create an produce
#           plots of a simulation of the PDE model that inject solute at
#           various different times and speeds.


import numpy as NP

import sys
sys.path.append('../../Utilities/')

from DistributionDataProcessor import DistributionDataProcessor;
from DistributionDataVisualiser import DistributionDataVisualiser;
from PDEModelSimulator import PDEModelSimulator;

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

end_time = 20;


#This function specifies the starting distribution of our Nanoparticle model.
sigma  = 1/14;
n0 = lambda x: NP.sqrt(0.5/NP.pi) / sigma * NP.exp(-1/2 * ((x - 1)/sigma)**2);


#The number of different injection functions that we wish to use in our
#simualtions out of the 6 availiable functions that have been implemented so
#far. This number should be a value between 0 and 6 (inclusive).
number_of_functions = 6;

#The size of the fonts to be used for the axis and legend labels of the output
#graphs.
label_font_size = 13;

#The size of the fonts to be used for the title of the output graphs.
title_font_size = 17;

##@brief    This example corresponds to the adding of solute to the
#           solution instantaneously at times t=5000, t=10000 and t=15000.
# @param t  This parameter corresponds to the current non-dimensionalised
#           version of time in the simulation that we wish to calculate
#           current $c_{\infty}(0)$ jump term for.
# @return   This function must return a number of type double that is >0,
#           which corresponds to the increase in concentration 
#           $c_{\infty}(0)$ caused by externally adding more solute to the
#           system.
def injectionFunction1(t):
    if(t>5000):
        if(t>10000):
            if(t>15000):
                return 30;
            else:
                return 20;
        else:
            return 10;
    else:
        return 0;
    
    
##@brief    This example corresponds to the adding of solute to the
#           solution instantaneously at times t=5000, and t=15000.
# @param t  This parameter corresponds to the current non-dimensionalised
#           version of time in the simulation that we wish to calculate
#           current $c_{\infty}(0)$ jump term for.
# @return   This function must return a number of type double that is >0,
#           which corresponds to the increase in concentration 
#           $c_{\infty}(0)$ caused by externally adding more solute to the
#           system.
def injectionFunction2(t):
    if(t>5000):
        if(t>15000):
            return 30;
        else:
            return 15;
    else:
        return 0;

##@brief    This example corresponds to smoothly adding solute to the reaction
#           between time 1000 and 16000.
# @param t  This parameter corresponds to the current non-dimensionalised
#           version of time in the simulation that we wish to calculate
#           current $c_{\infty}(0)$ jump term for.
# @return   This function must return a number of type double that is >0,
#           which corresponds to the increase in concentration 
#           $c_{\infty}(0)$ caused by externally adding more solute to the
#           system.
def injectionFunction3(t):
    if(t<1000):
        return 0;
    elif(t<16000):
        return 400*(t-1000)/15000
    else:
        return 400;

    
##@brief    This example corresponds to smoothly adding solute to the reaction
#           between time 100 and 1600.
# @param t  This parameter corresponds to the current non-dimensionalised
#           version of time in the simulation that we wish to calculate
#           current $c_{\infty}(0)$ jump term for.
# @return   This function must return a number of type double that is >0,
#           which corresponds to the increase in concentration 
#           $c_{\infty}(0)$ caused by externally adding more solute to the
#           system.
def injectionFunction4(t):
    if(t<100):
        return 0;
    elif(t<1600):
        return 600*(t-100)/1500
    else:
        return 600;


##@brief    This example corresponds to the gradual addition of solute 
#           to the solution between times 500 and 1000.
# @param t  This parameter corresponds to the current non-dimensionalised
#           version of time in the simulation that we wish to calculate
#           current $c_{\infty}(0)$ jump term for.
# @return   This function must return a number of type double that is >0,
#           which corresponds to the increase in concentration 
#           $c_{\infty}(0)$ caused by externally adding more solute to the
#           system.
def injectionFunction5(t):
    if(t<500):
        return 0;
    elif(t<1000):
        return 800*(t-500)/500
    else:
        return 800;

##@brief    This example corresponds to the gradual addition of solute to the
#           solution between times 1 and 10
# @param t  This parameter corresponds to the current non-dimensionalised
#           version of time in the simulation that we wish to calculate
#           current $c_{\infty}(0)$ jump term for.
# @return   This function must return a number of type double that is >0,
#           which corresponds to the increase in concentration 
#           $c_{\infty}(0)$ caused by externally adding more solute to the
#           system.
def injectionFunction6(t):
    if(t<1):
        return 0;
    elif(t<10):
        return 1000*(t-1)/9;
    else:
        return 1000;

#This array specifies the description strings in the legend for each line on
#the output graphs.
line_key_strings = ["3 sharp injections",
                    "2 sharp injections",
                    "smooth injection 1",
                    "smooth injection 2",
                    "smooth injection 3",
                    "smooth injection 4"];

#This array specifies the colour of each line on the output graphs.
line_colours = ["red",
                "green",
                "blue",
                "yellow",
                "purple",
                "orange"];

#This array specifies the type of each line on the output graphs.
line_styles = ["--",
               ":",
               "-",
               "-.",
               "-",
               "-"];

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
                        
#Retrieve the r and t values that will be used in the mesh grid for the
#simulation.
time_values = simulation_environment.getTValues()/3600;
r_values = simulation_environment.getRValues();


#Put the data into a processor so that we can create graphs an calculate
#some useful statistics.
data_processor = DistributionDataProcessor(r_values, time_values, 
                                           number_of_functions); 
for j in range(number_of_functions):
    data = simulation_environment.simulateModel(n0, injection_functions[j]);
    data_processor.addData(data, j);

#Create a visualiser for the simulation data.
data_visualiser = DistributionDataVisualiser("r (in m)", "N(r,t)", 
                                             line_key_strings[0:number_of_functions],
                                             line_colours[0:number_of_functions],
                                             line_styles[0:number_of_functions],
                                             label_font_size,
                                             title_font_size,
                                             data_processor);   

#Produce a graph of how the mean, variance and 3rd moment evolve over time.
data_visualiser.exportMeanTimeGraph(0, "Mean Plot Title", "t (in hours)",
                                    "Average r", "./mean_plot.png", True, 
                                    True, 0, 0);
data_visualiser.exportVarianceTimeGraph(1, "Variance Plot Title", 
                                        "t (in hours)", "Variance of r",
                                        "./variance_plot.png", True, True, 0,
                                        0);
data_visualiser.exportMomentTimeGraph(3, 2, "3rd Moment Plot Title",
                                      "t (in hours)", "Moment of r",
                                      "./moment_plot.png", True,
                                      True, 0, 0);
   
                                    
#Produce plots of every distribution at all time steps in the simulation.
data_visualiser.exportDistributionGraphAll("./Example_Distribution_Output_Images",
                                           False, True, 0, 0, True, 0, 0);