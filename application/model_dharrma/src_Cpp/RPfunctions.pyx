# cython: language_level=3
# distutils: language=c++
# cython: c_string_type=unicode, c_string_encoding=utf8

# Import necessary Cython declarations
from libcpp.vector cimport vector
from libcpp.string cimport string


# Import the C++ function declaration
cdef extern from "RPfunctions_src.cpp":

    cdef struct effFluid_Result:
       vector[double] kf, rhof, rhob

    cdef struct biotGassmann_Result:
        vector[double] Vp, Vs

    cdef struct hertzMindlin_Result:
        vector[double] KHM, muHM
    
    cdef struct hertzMindlin_trans_Result:
        vector[double] KHM, muHM

    cdef struct hillsAverage_Result:
        vector[double] mus, ks, rhos, nus

    effFluid_Result effFluid_src(vector[double] Sw, double kw, double ka, double rhow, double rhoa, vector[double] rhos, vector[string] soiltypes, vector[double] thicknesses, double dz)
    double fish_src(double vp, double vs)
    biotGassmann_Result biotGassmann_src(vector[double] KHM, vector[double] muHM, vector[double] ks, vector[double] kf, vector[double] rhob, vector[string] soiltypes, vector[double] thicknesses, double dz)
    hertzMindlin_Result hertzMindlin_src(vector[double] Swe, vector[double] z, vector[double] h, vector[double] rhob, double g, double rhoa, double rhow, vector[double] Ns, vector[double] mus, vector[double] nus, vector[double] fracs, int kk, vector[string] soiltypes, vector[double] thicknesses)
    hertzMindlin_trans_Result hertzMindlin_trans_src(vector[double] Swe, vector[double] z, vector[double] h, vector[double] rhob, vector[double] pression_profil, double g, double rhoa, double rhow, vector[double] Ns, vector[double] mus, vector[double] nus, vector[double] fracs, int kk, vector[string] soiltypes, vector[double] thicknesses)
    hillsAverage_Result hillsAverage_src(double mu_clay, double mu_silt, double mu_sand, double rho_clay, double rho_silt, double rho_sand, double k_clay, double k_silt, double k_sand, vector[string] soiltypes)



# Define Python wrapper functions
def effFluid(Sw, kw, ka, rhow, rhoa, rhos, soiltypes, thicknesses, dz):
    cdef effFluid_Result result = effFluid_src(Sw, kw, ka, rhow, rhoa, rhos, soiltypes, thicknesses, dz)
    return (result.kf, result.rhof, result.rhob)

def fish(vp, vs):
    cdef double sig = fish_src(vp, vs)
    return sig

def biotGassmann(KHM, muHM, ks, kf, rhob, soiltypes, thicknesses, dz):
    cdef biotGassmann_Result result = biotGassmann_src(KHM, muHM, ks, kf, rhob, soiltypes, thicknesses, dz)
    return (result.Vp, result.Vs)

def hertzMindlin(Swe, z, h, rhob, g, rhoa, rhow, Ns, mus, nus, fracs, kk, soiltypes, thicknesses):
    cdef hertzMindlin_Result result = hertzMindlin_src(Swe, z, h, rhob, g, rhoa, rhow, Ns, mus, nus, fracs, kk, soiltypes, thicknesses)
    return (result.KHM, result.muHM)

def hertzMindlin_trans(Swe, z, h, rhob, pression_profil, g, rhoa, rhow, Ns, mus, nus, fracs, kk, soiltypes, thicknesses):
    cdef hertzMindlin_trans_Result result = hertzMindlin_trans_src(Swe, z, h, rhob, pression_profil, g, rhoa, rhow, Ns, mus, nus, fracs, kk, soiltypes, thicknesses)
    return (result.KHM, result.muHM)

def hillsAverage(mu_clay, mu_silt, mu_sand, rho_clay, rho_silt, rho_sand, k_clay, k_silt, k_sand, soiltypes):
    cdef hillsAverage_Result result = hillsAverage_src(mu_clay, mu_silt, mu_sand, rho_clay, rho_silt, rho_sand, k_clay, k_silt, k_sand, soiltypes)
    return (result.mus, result.ks, result.rhos, result.nus)