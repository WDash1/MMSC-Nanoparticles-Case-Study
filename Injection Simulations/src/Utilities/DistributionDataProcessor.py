##  @brief  This file defines the DistributionDataProcessor class.


import numpy as NP
from StatisticsCalculator import StatisticsCalculator

##  @brief      This class is designed to provide the ability to manipulate and
#               axtract useful information (e.g. mean and variance) from the 
#               data produced by a series of simulations of the nanoparticle
#               PDE distribution model.
class DistributionDataProcessor:
    
    ##  @brief              The constructor for this class.
    #   @param r_values     This parameter should be an array of type double
    #                       that contains all the radius values used in the
    #                       mesh grid for simulations of the PDE model.    
    #   @param t_values     This parameter should be an array of type double
    #                       that contains the time value for each of the time
    #                       steps used in simulations of the PDE model.
    #   @param data_amount  This parameter should be an integer which is >0
    #                       and corresponds to the number of different
    #                       different numerical simulations of our model that 
    #                       we wish to process in this object. (e.g. if we
    #                       wish to compare two data sets one with solute
    #                       injection and one without, then this value should
    #                       be set to 2).
    def __init__(self, r_values, t_values, data_amount):
        self.__r_values = r_values;
        self.__t_values = t_values;
        self.__statistics_calculator = StatisticsCalculator(self.__r_values);
        self.__data_amount = data_amount;
        self.__raw_data = NP.zeros((self.__data_amount, len(self.__r_values), 
                                     len(self.__t_values))); 

    ##  @brief          This method may be used to add a set of results from
    #                   a numerical simulation of the PDE nanoparticle model.
    #   @param  N_data  This parameter must be an array of type double and size
    #                   nxm, where n is equal to the length of the r_values
    #                   parameter given to the constructor of this class and m
    #                   is equal to the length of the t_values parameter given
    #                   to the constructor of this class. Further the [i,j]th
    #                   entry of this array should correspond to the value of
    #                   the N distribution for a particle of size r_values[i]
    #                   at time t_values[j].
    #   @param  index   This parameter must be an integer which is >=0 and is
    #                   strictly less than the data_amount parameter given to
    #                   the constructor of this class. This parameter
    #                   corresponds to the index which we wish to use to
    #                   referance the data given to this method. It is
    #                   important to note that if the same index is used twice
    #                   then the most recent call to this method will overwrite
    #                   the data previously added using this index.        
    #   @note   In order to remove all the data previously added to this object
    #           by invoking this method, please use the clearData method.
    def addData(self, N_data, index):
        self.__raw_data[index]=N_data;

    ##  @brief      This method removes all the data added to this object by
    #               using the addData method of this class.
    def clearData(self):
        self.__raw_data = NP.zeros((self.data_amount, len(self.r_values), 
                                     len(self.t_values))); 

    ##  @brief          This method calculates the average value of a specified
    #                   distribution given to the constructor at each time step
    #                   in the simulation and returns the resulting values.
    #   @param index    This parameter must be an integer which is >=0 and is
    #                   strictly less than the data_amount parameter given to
    #                   the constructor of this class. This parameter
    #                   corresponds to the index (used when adding data to
    #                   this object via the addData method) of the
    #                   distribution data whose mean we wish to calculate.
    #   @return         This method returns an array of type double of the same
    #                   size as the t_values arrray given to the constructor
    #                   of this class. Further, the ith entry of the
    #                   aforementioned array corresponds to the average of the
    #                   simulation data stored at index "index" at time step
    #                   t_values[i].
    def calculateMeanValues(self, index):
        mean_function = lambda p, x: self.__statistics_calculator.calculateMean(p);
        mean_values = self.applyFunction(mean_function, index);
        return mean_values;
    
    
    ##  @brief      This method calculates the average of each of the 
    #               distributions encapsulated by this class at each time
    #               step in the simulation and returns the resulting values.
    #   @return     This method returns a two dimensional array of size nxm,
    #               where n is equal to the data_amount parameter given to the
    #               constructor of this class and m is equal to the size of the
    #               t_values array given to the constructor of this class.
    #               Further, the [i,j]th entry of the returned array
    #               corrsponds to the average of the ith distribution (in
    #               accordance with the addData method) at time step
    #               t_values[j].
    def calculateMeanValuesAll(self):
        mean_function = lambda p, x: self.__statistics_calculator.calculateMean(p);
        mean_values = self.applyFunctionToAll(mean_function);
        return mean_values;

    ##  @brief          This method calculates the variance of a specified
    #                   distribution given to the constructor at each time step
    #                   in the simulation and returns the resulting values.
    #   @param index    This parameter must be an integer which is >=0 and is
    #                   strictly less than the data_amount parameter given to
    #                   the constructor of this class. This parameter
    #                   corresponds to the index (used when adding data to
    #                   this object via the addData method) of the
    #                   distribution data whose variance we wish to calculate.
    #   @return         This method returns an array of type double of the same
    #                   size as the t_values arrray given to the constructor
    #                   of this class. Further, the ith entry of the
    #                   aforementioned array corresponds to the variance of the
    #                   simulation data stored at index "index" at time step
    #                   t_values[i].
    def calculateVarianceValues(self, index):
        variance_function = lambda p, x: self.__statistics_calculator.calculateVariance(p);
        variance_values = self.applyFunction(variance_function, index);
        return variance_values;


    
    ##  @brief      This method calculates the variance of each of the 
    #               distributions encapsulated by this class at each time
    #               step in the simulation and returns the resulting values.
    #   @return     This method returns a two dimensional array of size nxm,
    #               where n is equal to the data_amount parameter given to the
    #               constructor of this class and m is equal to the size of the
    #               t_values array given to the constructor of this class.
    #               Further, the [i,j]th entry of the returned array
    #               corrsponds to the variance of the ith distribution (in
    #               accordance with the addData method) at time step
    #               t_values[j].
    def calculateVarianceValuesAll(self):
        variance_function = lambda p, x: self.__statistics_calculator.calculateVariance(p);
        variance_values = self.applyFunctionToAll(variance_function);
        return variance_values;



    ##  @brief                  This method calculates the a designated moment
    #                           of a specified distribution given to the
    #                           constructor, at each time stepcin the
    #                           simulation and returns the resulting values.
    #   @param moment_number    This parameter must be an integer which is >0
    #                           and dictates the moment number that we wish to
    #                           calculate.
    #   @param index            This parameter must be an integer which is >=0
    #                           and is strictly less than the data_amount 
    #                           parameter given to the constructor of this
    #                           class. This parameter corresponds to the index
    #                           (used when adding data to this object via the
    #                           addData method) of the distribution data whose 
    #                           moment we wish to calculate.
    #   @return                 This method returns an array of type double of
    #                           the same size as the t_values arrray given to
    #                           the constructor of this class. Further, the ith
    #                           entry of the aforementioned array corresponds
    #                           to the variance of the simulation data stored 
    #                           at index "index" at time step t_values[i].
    def calculateMomentValues(self, moment_number, index):
        moment_function = lambda p, x: self.__statistics_calculator.calculateMoment(p, moment_number);
        moment_values = self.applyFunction(moment_function, index);
        return moment_values;


    
    ##  @brief                  This method calculates the a specified moment
    #                           of each of the distributions encapsulated by 
    #                           this class at each time step in the simulation
    #                           and returns the resulting values.
    #   @param moment_number    This parameter must be an integer which is >0
    #                           and corresponds to the moment number of
    #                           the distributions that we wish to compute.
    #   @return                 This method returns a two dimensional array of
    #                           size nxm, where n is equal to the data_amount
    #                           parameter given to the constructor of this
    #                           class and m is equal to the size of the 
    #                           t_values array given to the constructor of 
    #                           this class. Further, the [i,j]th entry of the
    #                           returned array corrsponds to the specified
    #                           moment of the ith distribution (in accordance
    #                           with the addData method) at time step
    #                           t_values[j].
    def calculateMomentValuesAll(self, moment_number):
        moment_function = lambda p, x: self.__statistics_calculator.calculateMoment(p, moment_number);
        moment_values = self.applyFunctionToAll(moment_function);
        return moment_values;

    
    
    ##  @brief                  This method applies a given function
    #                           to a specified distribution given to the
    #                           constructor, at each time step in the
    #                           simulation and returns the resulting values.
    #   @param func             This parameter must be a function handle which
    #                           takes two arrays with the same length as the
    #                           r_values array given to the constructor and
    #                           returns a value of type double. (The first
    #                           argument given to this function is the
    #                           value of the distirbution N at a particular
    #                           time step and the second is the r_values array
    #                           that dictates the size of the nanoparticle
    #                           for each entry in the N distribution array).
    #   @param index            This parameter must be an integer which is >=0
    #                           and is strictly less than the data_amount 
    #                           parameter given to the constructor of this
    #                           class. This parameter corresponds to the index
    #                           (used when adding data to this object via the
    #                           addData method) of the distribution data to
    #                           which we wish to apply the function func at
    #                           each time step.
    #   @return                 This method returns an array of type double of
    #                           the same size as the t_values arrray given to
    #                           the constructor of this class. Further, the ith
    #                           entry of the aforementioned array corresponds
    #                           to applying the function "func" to the
    #                           simulation data stored at index "index" at time
    #                           step t_values[i].
    def applyFunction(self, func, index):
        output_values = NP.empty(len(self.__t_values));
        for j in range(len(self.__t_values)):
            output_values[j] = func(self.__raw_data[index][:,j], self.__r_values);
        return output_values;
    
    ##  @brief          This method applied a designated function to each of 
    #                   the distributions encapsulated by this class at each
    #                   time step in the simulation and returns the resulting
    #                   values.
    #   @param func     This parameter must be a function handle which
    #                   takes two arrays with the same length as the
    #                   r_values array given to the constructor and
    #                   returns a value of type double. (The first
    #                   argument given to this function is the
    #                   value of the distirbution N at a particular
    #                   time step and the second is the r_values array
    #                   that dictates the size of the nanoparticle
    #                   for each entry in the N distribution array).
    #   @return         This method returns a two dimensional array of size
    #                   nxm, where n is equal to the data_amount parameter
    #                   given to the constructor of this class and m is equal
    #                   to the size of the t_values array given to the
    #                   constructor of this class. Further, the [i,j]th entry 
    #                   of the returned array corresponds to the result of 
    #                   the function "func" applied to the ith distribution (in
    #                   accordance with the addData method) at time step
    #                   t_values[j].
    def applyFunctionToAll(self, func):
        output_values = NP.zeros((self.__data_amount, len(self.__t_values)));
        for index in range(self.__data_amount):
            output_values[index] = self.applyFunction(func, index);
        return output_values;
    
    
    ##  @brief              This method may be used to access a paticular
    #                       distribution encapsulated by this class at a
    #                       designated time step in the simulation.
    #   @param index        This parameter must be an integer which is >=0 and
    #                       is strictly less than the data_amount parameter
    #                       given to the constructor of this class. This
    #                       parameter corresponds to the index of the data set
    #                       whose distribution values we wish to access (note:
    #                       this index must be consistent with the index we
    #                       used for adding data to this class by using the
    #                       addData method of this class).
    #   @param time_index   This prameter must be an integer which is >=0 and
    #                       is strictly less than the t_values array given to
    #                       the constructor of this class. This parameter 
    #                       dictates the index of the entry in the t_values
    #                       array for which we wish to retrieve the values of
    #                       the specified distribtuion.
    #   @return             This method returns an array of the same size as
    #                       the r_values array given to the constructor of
    #                       this class. Further, the ith entry of this array
    #                       corresponds to the value of our distribution N for
    #                       r value: r_values[i] at time: t_values[time_index].
    def getDistribution(self, index, time_index):
        return NP.array(self.__raw_data)[index][:,time_index];


    ##  @brief              This method may be used to access each of the
    #                       distributions encapsulated by this class at 
    #                       a specified time step in the simulation.
    #   @param time_index   This prameter must be an integer which is >=0 and
    #                       is strictly less than the t_values array given to
    #                       the constructor of this class. This parameter 
    #                       dictates the index of the entry in the t_values
    #                       array for which we wish to retrieve the values of
    #                       the distributions encapsulated by this class.
    #   @return             This method returns a two dimensional array of size
    #                       nxm, where n is equal to the data_amount parameter
    #                       given to the constructor of this class and m is
    #                       equal to the length of the r_values array given to
    #                       the constructor of this class. Further, the [i,j]th
    #                       entry of this array corresponds to the value of the
    #                       distribution with index i (in accordance with the
    #                       addData method) for nanoparticle radius value
    #                       r_values[j].
    def getDistributionAll(self, time_index):
        output = list(map(lambda i: self.__raw_data[i][:,time_index], 
                          list(range(self.__data_amount))));
        return output;
    
    
    
    ##  @brief      This method is used to access the r_values array given to
    #               the constructor of this class.
    #   @return     This method returns the r_values array given to the
    #               constructor of this class.    
    def getRValues(self):
        return self.__r_values;
        
    ##  @brief      This method is used to access the t_values array given to
    #               the constructor of this class.
    #   @return     This method returns the t_values array given to the 
    #               constructor of this class.
    def getTValues(self):
        return self.__t_values;