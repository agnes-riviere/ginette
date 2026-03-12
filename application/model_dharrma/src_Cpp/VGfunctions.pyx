# cython: language_level=3
# distutils: language=c++
# cython: c_string_type=unicode, c_string_encoding=utf8

# Import necessary Cython declarations
from libcpp.vector cimport vector
from libcpp.string cimport string



# Import the C++ function declaration
cdef extern from "VGfunctions_src.cpp":

    struct vanGen_Result:
        vector[double] h, Sw, Swe

    struct selectSoilType_Result:
        double wsand, wclay, wsilt, phi, alpha, nvg, theta, Swr

    vanGen_Result vanGen_src(vector[double] z, double WT, vector[string] soiltypes, vector[double] thicknesses)
    selectSoilType_Result selectSoilType_src(string soiltype)



# Define Python wrapper functions
def vanGen(vector[double] z, double WT, vector[string] soiltypes, vector[double] thicknesses):
    cdef vanGen_Result result = vanGen_src(z, WT, soiltypes, thicknesses)
    return result.h, result.Sw, result.Swe

def selectSoilType(string soiltype):
    cdef selectSoilType_Result result = selectSoilType_src(soiltype)
    return result.wsand, result.wclay, result.wsilt, result.phi, result.alpha, result.nvg, result.theta, result.Swr 
