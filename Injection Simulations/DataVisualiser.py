import pylab
import numpy

## @brief   This class provides the ability to  visualise and export data
#           which has been produced by numerical simulations of our ODE
#           and PDE equations.
class DataVisualiser:
   
    ## @brief   The constructor for this class.
    #  @param figure_num    This parameter should be an integer which is >0
    #                       that corresponds to the figure number associated
    #                       with graphs produced by this class.
    #  @param x_lower_lim   This parameter should be a number of type double
    #                       that is strictly less than the x_upper_lim
    #                       paramter and corresponds to the smallest value
    #                       displayed on the x-axis of graphs produced by this
    #                       class.
    #  @param x_upper_lim   This parameter should be a number of type double
    #                       that is strictly greater than the x_lower_lim
    #                       parameter and corresponds to the largest value
    #                       displayed on the x-axis of the graphs produced by
    #                       this class.
    #  @param x_label       This parameter should be a string that dictates the
    #                       label that will be displayed on the horizontal
    #                       axis of the graphs produced by this class.
    #  @param y_lower_lim   This parameter should be a number of type double
    #                       that is strictly less than the y_upper_lim
    #                       parameter and corresponds to the smallest value
    #                       displayed on the y-axis of the graphs produced by
    #                       this class.
    #  @param y_upper_lim   This parameter should be a number of type double
    #                       that is strictly greater than the y_lower_lim
    #                       parameter and corresponds to the largest value
    #                       displayed on the y-axis of the graphs produced by
    #                       this class.
    #  @param y_label       This parameter should be a string that dictates the
    #                       label that will be displayed on the vertical
    #                       axis of the graphs produced by this class.
    #  @param key_strings   This parameter should be an array of strings that
    #                       dictate the discriptions shown for each data set,
    #                       in the legend of the graphs generated by this 
    #                       class. It is important to note that this array
    #                       must be the same length as the "line_colours"
    #                       array.
    #  @param line_colours  This parameter should be an array of strings that
    #                       dictate the colours to be used for each data set
    #                       displayed in the graphs generated by this class
    #                       (e.g. this parameter may be set to ["red",
    #                       "green", "blue", "yellow", "purple", "orange"].
    #                       It is important to note that this array must be
    #                       the same length as the "key_strings" array.
    def __init__(self, figure_num, 
                 x_lower_lim, x_upper_lim, x_label, 
                 y_lower_lim, y_upper_lim, y_label, 
                 key_strings, line_colours):

        self.key_strings = key_strings;
        self.line_colours = line_colours;
        self.figure_num = figure_num;
        
        self.x_lower_lim = x_lower_lim;
        self.x_upper_lim = x_upper_lim;
        
        self.y_lower_lim = y_lower_lim;
        self.y_upper_lim = y_upper_lim;
        
        self.x_label = x_label;
        self.y_label = y_label;
        
        self.data_set_amount = len(key_strings); 
        
        self.plot_data = numpy.empty(self.data_set_amount, dtype=object);
        
        
        
    ##  @brief              This method of the class may be used to add a set
    #                       of data to be displayed on the graph produced by 
    #                       this graph.
    #   @param x_values     This parameter must be an array of type double that
    #                       corresponds to the x-coordinate of each of the 
    #                       points that we wish to display in the graph
    #                       produced by this class. It is important to note
    #                       that this array must be the same length as the 
    #                       "y_values" parameter.
    #   @param y_values     This parameter must be an array of type double that
    #                       corresponds to the y-coordinate of each of the
    #                       points that we wish to display in the graph
    #                       produced by this class. It is important to note
    #                       that this array must be the same length as the
    #                       "x_values" parameter.
    #   @param index        This parameter must be an integer that is >=0 that
    #                       dictates which colour and description settings
    #                       given to the construct should be associated with
    #                       the x and y data given to this method.
    #                       It is important to note that this index value must
    #                       be strictly less than the length of the arrays
    #                       "key_strings" and "line_colours" given to the
    #                       constructor of this class.
    #
    def addData(self, x_values, y_values, index):
            self.plot_data[index] = [x_values, y_values];
   
    # @brief                    This method produces a line graph from the data
    #                           which has been added to this class using the
    #                           settings given to the constructor. Further, 
    #                           this method will export the aforementioned
    #                           graph to a specified image file.
    # @param current_title      This parameter is a string which dictates the
    #                           title that will be displayed on the graph
    #                           generated by this method.
    # @param output_fileapth    This parameter specifies the path of the image
    #                           file to which we wish to write the generated
    #                           graph image file to.    
    # @param use_log_scale      This parameter should be a boolean that
    #                           dictates whether the provided data is
    #                           displayed on a logarithmic scale in the 
    #                           resulting graph.
    def exportGraph(self, current_title, output_filepath, use_log_scale):
        pylab.figure(self.figure_num);
        pylab.title(current_title);
        pylab.xlabel(self.x_label);
        pylab.ylabel(self.y_label);
        pylab.xlim(self.x_lower_lim, self.x_upper_lim);
        pylab.ylim(self.y_lower_lim, self.y_upper_lim);
        
        for i in range (len(self.line_colours)):
            if(use_log_scale==True):
                pylab.semilogy(self.plot_data[i][0], self.plot_data[i][1], 'k--', 
                    label = self.key_strings[i], color = self.line_colours[i]);
            else:
                pylab.plot(self.plot_data[i][0], self.plot_data[i][1], 'k--', 
                    label = self.key_strings[i], color = self.line_colours[i]);
        
        pylab.legend(fontsize = 10);
        pylab.savefig(output_filepath, bbox_inches = 'tight');
        pylab.close(self.figure_num);
    
    
    # @brief    This method will erase all data previosly added to this class.
    def clearData(self):
        self.plot_data = numpy.empty(self.data_set_amount, dtype=object);