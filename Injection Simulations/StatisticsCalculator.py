import numpy as NP

## @brief   This class provides a method for calculating the mean, variance
#           and arbitrary moments of a given distribution function.
class StatisticsCalculator:
    ## @brief           The constructor for this class.
    #  @param x_values  This parameter should be an array of doubles that
    #                   corresponds to the x values associated with each of
    #                   the probability weights given to the
    #                   calculateMoment method of this class.
    def __init__(self, x_values):
        self.__x_values = x_values;

    ## @brief           This method normalises a given distribution function so
    #                   that its entries all sum to 1.
    #  @param p_values  This parameter should be an array of doubles which are
    #                   all >=0 and whose size is equal
    #                   to that of the "x_values" array given to the
    #                   constructor of this class. This array corresponds to
    #                   the probability values for our distribution that we
    #                   wish to normalise.
    #  @return          This method returns a numpy array, which contains
    #                   scaled versions of each of the entries in the
    #                   "p_values" parameter, so that all the entries in the
    #                   vector sum to 1.
    def noramliseData(self, p_values):
        return NP.array(p_values)*(1/sum(p_values));

    ## @brief                   This method calculates a given moment
    #                           associated with a specified list of
    #                           probability values.
    #  @param p_values          This parameter must be an array of type double,
    #                           which is the same length as the x_values 
    #                           array given to the constructor of this class. 
    #                           The entries in this array should all be
    #                           >=0 and correspond to the probability of each
    #                           x value in the distribution occuring. 
    #                           It is important to note that if the entries in
    #                           this array to not sum to 1, then this method
    #                           will normalise the array by using the
    #                           normaliseData method of this class.
    #  @param moment_number     This parameter must be an integer which is
    #                           strictly greater than zero, that dictates
    #                           the moment number of the provided probability
    #                           distribution that we wish to compute.
    #  @return                  This method will return a value of type double
    #                           that corresponds to the designated moment of
    #                           the provided probability distribution.
    def calculateMoment(self, p_values, moment_number):
        distribution_data = NP.array(p_values)*(1/sum(p_values));
        moment_value = NP.dot(distribution_data, NP.power(self.__x_values, 
                                                          moment_number));
        return moment_value;

    ## @brief                   This method calculates the mean of the 
    #                           distribution encapsulated by this class,
    #                           for a specified list of probability values.
    #  @param p_values          This parameter must be an array of type double,
    #                           which is the same length as the x_values 
    #                           array given to the constructor of this class. 
    #                           The entries in this array should all be
    #                           >=0 and correspond to the probability of each
    #                           x value in the distribution occuring. 
    #                           It is important to note that if the entries in
    #                           this array to not sum to 1, then this method
    #                           will normalise the array by using the
    #                           normaliseData method of this class.
    # @return                   This method will return a value of type double
    #                           that corresponds to the mean of
    #                           the provided probability distribution.            
    def calculateMean(self, p_values):
        return self.calculateMoment(p_values, 1);
    
    
    ## @brief                   This method calculates the variance of the
    #                           distribution encapsulated by this class,
    #                           for a specific list of probability values.
    #  @param p_values          This parameter must be an array of type double,
    #                           which is the same length as the x_values 
    #                           array given to the constructor of this class. 
    #                           The entries in this array should all be
    #                           >=0 and correspond to the probability of each
    #                           x value in the distribution occuring. 
    #                           It is important to note that if the entries in
    #                           this array to not sum to 1, then this method
    #                           will normalise the array by using the
    #                           normaliseData method of this class.
    #  @return                  This method will return a value of type double
    #                           that corresponds to the variance of the
    #                           the provided probability distribution.
    def calculateVariance(self, p_values):
        return self.calculateMoment(p_values, 2) - self.calculateMean(p_values)**2;