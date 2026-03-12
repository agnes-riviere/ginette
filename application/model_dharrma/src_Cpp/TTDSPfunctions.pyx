# cython: language_level=3
# distutils: language=c++
# cython: c_string_type=unicode, c_string_encoding=utf8

# Import necessary Cython declarations
from libcpp.string cimport string
from libcpp.vector cimport vector
import numpy as np



# Import the C++ function declaration
cdef extern from "TTDSPfunctions_src.cpp":
    string writeVelocityModel_src(vector[double] thk, vector[double] vp, vector[double] vs,  vector[double] rho, string substratum, int n_layers_substratum)
    vector[double] firstArrival_src(vector[double] thk, vector[double] vv, vector[double] Xdata, double trig)
    

# Define Python wrapper functions
def writeVelocityModel(thk, vp, vs, rho, substratum, n_layers_substratum):
    cdef string velocity_model = writeVelocityModel_src(thk, vp, vs, rho, substratum, n_layers_substratum)
    return velocity_model

def firstArrival(thk, vv, Xdata, trig):
    cdef vector[double] thk_cpp = thk
    cdef vector[double] vv_cpp = vv
    cdef vector[double] Xdata_cpp = Xdata
    cdef vector[double] result = firstArrival_src(thk_cpp, vv_cpp, Xdata_cpp, trig)
    Thod = np.array(result)
    return Thod


# Define a Python function
def readDispersion(string gpdc_output_string):
    '''
    Reads the dispersion data from the output string of the GPDC program.

    Parameters:
    gpdc_output_string (str): The output string of the GPDC program

    Outputs:
    dispersion_data (list): A list of numpy arrays, each containing the dispersion data for a single mode
    n_modes (int): The number of computed modes in the dispersion data
    '''
    # cleaned output string without comments and empty lines
    gpdc_output_string = '\n'.join(line for line in gpdc_output_string.split('\n') if not line.strip().startswith('#') and line.strip())
    
    # array with all dispersion data not separated by mode
    gpdc_output_array = np.fromstring(gpdc_output_string, dtype=float, sep=' ').reshape(-1, 2)

    # Convert dispersion data to a list of numpy arrays for each mode
    dispersion_data = []
    current_array = []
    for freq, vel in gpdc_output_array:
        if len(current_array) == 0 or freq >= current_array[-1][0]:
            current_array.append((freq, 1.0/vel))
        else:
            dispersion_data.append(np.array(current_array))
            current_array = [(freq, 1.0/vel)]
    if current_array:
        dispersion_data.append(np.array(current_array))
    return dispersion_data, len(dispersion_data)