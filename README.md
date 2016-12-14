# serpentine-model

Model written in R simulates a Serpentine landscape occupied by plants competing for space; Project for fullfillment of requirements of a Bachelors of Arts, Bennington College. 

Thesis paper which documents this model can be found at: 

Running the model:

parameter_graphing is the code which parameterizes and runs multiple simulaltions.
The code is dependent on the Serpentine_model file, whose location is referenced in the  parameter_graphing code and thus needs to be specified before attempting to run the code. 

Functions for the code must also be run. These functions are located in the  miekes.model_functions directory. Each function is stored individually in a text file. Some functions are not used in the current model (these should be trimmed from the project or archived).

You will also need to specifiy where data is to be stored. This should be specified in the parameter_graphing code.

data_analysis.Rmd documents and explains the code used for data analysis done with the model results. This is also available as a PDF. 

fract_trial_3_6spec.txt contains the simulation data used in the data analysis.
