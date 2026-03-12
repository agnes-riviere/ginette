#include <cmath>
#include <vector>
#include "VGfunctions_src.h"
#include <iostream>



struct effFluid_Result {
    std::vector<double> kf, rhof, rhob;
};

struct biotGassmann_Result {
    std::vector<double> Vp, Vs;
};

struct hertzMindlin_Result {
    std::vector<double> KHM, muHM;
};

struct hertzMindlin_trans_Result {
    std::vector<double> KHM, muHM;
};

struct hillsAverage_Result {
    std::vector<double> mus, ks, rhos, nus;
};



hillsAverage_Result hillsAverage_src(const double& mu_clay,
                                     const double& mu_silt,
                                     const double& mu_sand,
                                     const double& rho_clay,
                                     const double& rho_silt,
                                     const double& rho_sand,
                                     const double& k_clay,
                                     const double& k_silt,
                                     const double& k_sand,
                                     const std::vector<std::string>& soiltypes) {
    // ==========================================
    // Computes the effective properties of the solid grains from its constituents
    //
    // Parameters:
    // double mu_clay: Shear modulus of clay [GPa]
    // double mu_silt: Shear modulus of silt [GPa]
    // double mu_sand: Shear modulus of sand [GPa]
    // double rho_clay: Density of clay [Kg/m3]
    // double rho_silt: Density of silt [Kg/m3]
    // double rho_sand: Density of sand [Kg/m3]
    // double k_clay: Bulk modulus of clay [GPa]
    // double k_silt: Bulk modulus of silt [GPa]
    // double k_sand: Bulk modulus of sand [GPa]
    // vector[string] soiltypes: soil types of each layer in the profile
    //
    // Outputs:
    // hillsAverage_Result result: structure with the following fields
    // vector [double] mus: Shear moduli of grain of each layer [Pa]
    // vector[double] ks: Bulk moduli of grain of each layer [Pa]
    // vector[double] rhos: Densities of grain of each layer [Kg/m3]
    // vector[double] nus: Poisson's ratios of each layer [-]
    //
    // Programmers:
    // Firstly implemented in Matlab in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // Adaptated for multi-layers by J. Cunha Teixeira in 2024/04
    // ==========================================

    std::vector<double> mus(soiltypes.size()); // Shear moduli of grain of each layer
    std::vector<double> ks(soiltypes.size()); // Bulk moduli of grain of each layer
    std::vector<double> rhos(soiltypes.size()); // Densities of grain of each layer
    std::vector<double> nus(soiltypes.size()); // Poisson's ratios of each layer

    for (size_t j = 0; j < soiltypes.size(); ++j) {
        std::string_view soilType = soiltypes[j]; // Soil type for current layer

        selectSoilType_Result soil = selectSoilType_src(soilType); // Soil properties for current layer

        // Shear moduli of grain [Pa]
        mus[j] = 0.5 * (1 / (soil.wclay / mu_clay + soil.wsilt / mu_silt + soil.wsand / mu_sand) + (soil.wclay * mu_clay + soil.wsilt * mu_silt + soil.wsand * mu_sand)) * 1e9;
        
        // Bulk moduli of grain [Pa]
        ks[j] = 0.5 * (1 / (soil.wclay / k_clay + soil.wsilt / k_silt + soil.wsand / k_sand) + (soil.wclay * k_clay + soil.wsilt * k_silt + soil.wsand * k_sand)) * 1e9;
        
        // Densities of grain [Kg/m3]
        rhos[j] = soil.wclay * rho_clay + soil.wsilt * rho_silt + soil.wsand * rho_sand;
        
        // Poisson's ratios
        nus[j] = (3 * ks[j] - 2 * mus[j]) / (2 * (3 * ks[j] + mus[j]));
    }

    hillsAverage_Result result; // Structure to store the results
    result.mus = mus;
    result.ks = ks;
    result.rhos = rhos;
    result.nus = nus;

    return result;
}



effFluid_Result effFluid_src(const std::vector<double>& Sws,
                             const double& kw,
                             const double& ka,
                             const double& rhow,
                             const double& rhoa,
                             const std::vector<double>& rhos,
                             const std::vector<std::string>& soiltypes,
                             const std::vector<double>& thicknesses,
                             const double& dz) {
    // ==========================================
    // Computes effective fluid properties and bulk density
    //
    // Parameters:
    // vector[double] Sws: total saturation over depth [-]
    // double kw: water bulk modulus [Pa]
    // double ka: air bulk modulus [Pa]
    // double rhow: water density [Kg/m3]
    // double rhoa: air density [Kg/m3]
    // vector[double] rhos: Densities of grain of each layer [Kg/m3]
    // vector[string] soiltypes: soil types of each layer in the profile
    // vector[double] thicknesses: thickness of each layer [m]
    // double dz: depth discretization [m]
    //
    // Outputs:
    // effFluid_Result result: structure with the following fields
    // vector [double] kf: effective compressibility over depth [Pa-1]
    // vector[double] rhof: effective fluid density over depth [Kg/m3]
    // vector[double] rhob: bulk density over depth [Kg/m3]
    //
    // Programmers:
    // Firstly implemented in Matlab in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // Adaptated for multi-layers by J. Cunha Teixeira in 2024/04
    // ==========================================

    size_t num_points = Sws.size(); // Number of cells in the profile

    std::vector<double> kf(num_points); // Effective compressibility over depth
    std::vector<double> rhof(num_points); // Effective fluid density over depth
    std::vector<double> rhob(num_points); // Bulk density over depth

    int start = 0; // Depth start index for fisrt layer

    for (size_t j = 0; j < soiltypes.size(); ++j) {
        std::string_view soilType = soiltypes[j]; // Soil type for current layer
        double thickness = thicknesses[j]; // Thickness of current layer

        selectSoilType_Result soil = selectSoilType_src(soilType); // Soil properties for current layer

        int end = round(thickness / dz) + start; // Depth end index for current layer

        for (int i = start; i < end; ++i) {
            kf[i] = 1.0 / (Sws[i] / kw + (1.0 - Sws[i]) / ka); // Effective compressibility (Woods formula)
            rhof[i] = Sws[i] * rhow + (1.0 - Sws[i]) * rhoa; // Effective fluid density
            rhob[i] = (1.0 - soil.phi) * rhos[j] + soil.phi * rhof[i]; // Bulk density
        }
        start = end; // Update start index for next layer
    }

    effFluid_Result result; // Structure to store the results
    result.kf = kf;
    result.rhof = rhof;
    result.rhob = rhob;

    return result;
}



hertzMindlin_Result hertzMindlin_src(const std::vector<double>& Swe,
                                     const std::vector<double>& z,
                                     const std::vector<double>& h,
                                     const std::vector<double>& rhob,
                                     const double& g,
                                     const double& rhoa,
                                     const double& rhow,
                                     const std::vector<double>& Ns,
                                     const std::vector<double>& mus,
                                     const std::vector<double>& nus,
                                     const std::vector<double>& fracs,
                                     const int& kk,
                                     const std::vector<std::string>& soiltypes,
                                     const std::vector<double>& thicknesses) {
    // ==========================================
    // Hertz-Mindlin model
    //
    // Parameters:
    // vector[double] Swe: effective wetting phase saturation over depth [-]
    // vector[double] z: depths of the profile [m]
    // vector[double] h: pressure head over depth [m]
    // vector[double] rhob: bulk density over depth [Kg/m3]
    // double g: gravity [m/s2]
    // double rhoa: air density [Kg/m3]
    // double rhow: water density [Kg/m3]
    // vector[double] Ns: coordination number of each layer [-]
    // vector[double] mus: Shear moduli of grain of each layer [GPa]
    // vector[double] nus: Poisson's ratios of each layer [-]
    // vector[double] fracs: fraction of non-slipping grains of each layer [-]
    // int kk: type of effective stress
    // vector[string] soiltypes: soil types of each layer in the profile
    // vector[double] thicknesses: thickness of each layer [m]
    //
    // Outputs:
    // hertzMindlin_Result result: structure with the following fields
    // vector [double] KHM: effective bulk modulus over depth [Pa]
    // vector[double] muHM: effective shear modulus over depth [Pa]
    //
    // Programmers:
    // Firstly implemented in Matlab in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // Adaptated for multi-layers by J. Cunha Teixeira in 2024/04
    // ==========================================

    size_t num_points = Swe.size(); // Number of cells in the profile

    std::vector<double> KHM(num_points); // Effective bulk modulus over depth
    std::vector<double> muHM(num_points); // Effective shear modulus over depth
    
    std::vector<double> Pe(num_points); // Overburden stress over depth

    double dz = fabs(z[1] - z[0]); // Depth discretization

    int start = 0; // Depth start index for fisrt layer

    double sigma_prev = 0; // Previous overburden stress

    for (size_t j = 0; j < soiltypes.size(); ++j) {
        std::string_view soilType = soiltypes[j]; // Soil type for current layer
        double thickness = thicknesses[j]; // Thickness of current layer

        selectSoilType_Result soil = selectSoilType_src(soilType); // Soil properties for current layer

        int end = round(thickness / dz) + start; // Depth end index for current layer

        for (int i = start; i < end; ++i) {

            double sigma = rhob[i] * g * dz + sigma_prev; // Overburden stress
            sigma_prev = sigma;

            double dep = fabs(z[i]); // Depth (positive value)
            double St = (h[i] <= 0) ? 0 : Swe[i] * rhow * g * h[i]; // Suction term | Null when fully saturated
            double Pa = rhoa * g * dep; // Air pressure
            double Pe; // Effective stress
            
            if (kk == 1) {
                Pe = 101325;
            } else if (kk == 2) {
                Pe = sigma - Pa;
            } else if (kk == 3) {
                Pe = sigma - Pa + St;
            } else {
                throw std::invalid_argument("Invalid value for kk");
            }

            KHM[i] = std::pow((Ns[j]*Ns[j] * std::pow((1 - soil.phi), 2) * mus[j]*mus[j] * Pe) / (18 * M_PI*M_PI * std::pow((1 - nus[j]), 2)), 1.0 / 3.0); // Effective bulk modulus
            muHM[i] = ((2 + 3 * fracs[j] - (1 + 3 * fracs[j]) * nus[j]) / (5 * (2 - nus[j]))) * std::pow((3 * Ns[j]*Ns[j] * std::pow((1 - soil.phi), 2) * mus[j]*mus[j] * Pe) / (2 * M_PI*M_PI * std::pow((1 - nus[j]), 2)), 1.0 / 3.0); // Effective shear modulus
        }
        start = end; // Update start index for next layer
    }

    hertzMindlin_Result result; // Structure to store the results
    result.KHM = KHM;
    result.muHM = muHM;
    
    return result;
}

hertzMindlin_trans_Result hertzMindlin_trans_src(const std::vector<double>& Swe,
                                     const std::vector<double>& z,
                                     const std::vector<double>& h,
                                     const std::vector<double>& rhob,
                                     const std::vector<double>& pression_profil,
                                     const double& g,
                                     const double& rhoa,
                                     const double& rhow,
                                     const std::vector<double>& Ns,
                                     const std::vector<double>& mus,
                                     const std::vector<double>& nus,
                                     const std::vector<double>& fracs,
                                     const int& kk,
                                     const std::vector<std::string>& soiltypes,
                                     const std::vector<double>& thicknesses) {
    // ==========================================
    // Hertz-Mindlin model
    //
    // Parameters:
    // vector[double] Swe: effective wetting phase saturation over depth [-]
    // vector[double] z: depths of the profile [m]
    // vector[double] h: pressure head over depth [m]
    // vector[double] rhob: bulk density over depth [Kg/m3]
    // double g: gravity [m/s2]
    // double rhoa: air density [Kg/m3]
    // double rhow: water density [Kg/m3]
    // vector[double] Ns: coordination number of each layer [-]
    // vector[double] mus: Shear moduli of grain of each layer [GPa]
    // vector[double] nus: Poisson's ratios of each layer [-]
    // vector[double] fracs: fraction of non-slipping grains of each layer [-]
    // int kk: type of effective stress
    // vector[string] soiltypes: soil types of each layer in the profile
    // vector[double] thicknesses: thickness of each layer [m]
    //
    // Outputs:
    // hertzMindlin_Result result: structure with the following fields
    // vector [double] KHM: effective bulk modulus over depth [Pa]
    // vector[double] muHM: effective shear modulus over depth [Pa]
    //
    // Programmers:
    // Firstly implemented in Matlab in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // Adaptated for multi-layers by J. Cunha Teixeira in 2024/04
    // ==========================================

    size_t num_points = Swe.size(); // Number of cells in the profile

    std::vector<double> KHM(num_points); // Effective bulk modulus over depth
    std::vector<double> muHM(num_points); // Effective shear modulus over depth
    
    std::vector<double> Pe(num_points); // Overburden stress over depth

    double dz = fabs(z[1] - z[0]); // Depth discretization

    int start = 0; // Depth start index for fisrt layer

    double sigma_prev = 0; // Previous overburden stress

    for (size_t j = 0; j < soiltypes.size(); ++j) {
        std::string_view soilType = soiltypes[j]; // Soil type for current layer
        double thickness = thicknesses[j]; // Thickness of current layer

        selectSoilType_Result soil = selectSoilType_src(soilType); // Soil properties for current layer

        int end = round(thickness / dz) + start; // Depth end index for current layer

        for (int i = start; i < end; ++i) {

            double sigma = rhob[i] * g * dz + sigma_prev; // Overburden stress
            sigma_prev = sigma;
            double St;

            double dep = fabs(z[i]); // Depth (positive value)
            if (kk == 4) {
                St = (h[i] <= 0) ? 0 : Swe[i] * abs(pression_profil[i]); // Suction term | Null when fully saturated Utilisation de la pression ginette
                // std::cout << "Bonjour le monde !" << std::endl;
            } else {
                St = (h[i] <= 0) ? 0 : Swe[i] * rhow * g * h[i]; // Suction term | Null when fully saturated
            }

            double Pa = rhoa * g * dep; // Air pressure 
            double Pe; // Effective stress (pas une pression) sigma
            
            if (kk == 1) {
                Pe = 101325;
            } else if (kk == 2) {
                Pe = sigma - Pa;
            } else if (kk == 3) {
                Pe = sigma - Pa + St;
            } else if (kk == 4) {
                Pe = sigma - Pa + St;
            } else {
                throw std::invalid_argument("Invalid value for kk");
            }

            KHM[i] = std::pow((Ns[j]*Ns[j] * std::pow((1 - soil.phi), 2) * mus[j]*mus[j] * Pe) / (18 * M_PI*M_PI * std::pow((1 - nus[j]), 2)), 1.0 / 3.0); // Effective bulk modulus
            muHM[i] = ((2 + 3 * fracs[j] - (1 + 3 * fracs[j]) * nus[j]) / (5 * (2 - nus[j]))) * std::pow((3 * Ns[j]*Ns[j] * std::pow((1 - soil.phi), 2) * mus[j]*mus[j] * Pe) / (2 * M_PI*M_PI * std::pow((1 - nus[j]), 2)), 1.0 / 3.0); // Effective shear modulus
        }
        start = end; // Update start index for next layer
    }

    hertzMindlin_trans_Result result; // Structure to store the results
    result.KHM = KHM;
    result.muHM = muHM;
    
    return result;
}



biotGassmann_Result biotGassmann_src(const std::vector<double>& KHM,
                                     const std::vector<double>& muHM,
                                     const std::vector<double>& ks,
                                     const std::vector<double>& kf,
                                     const std::vector<double>& rhob,
                                     const std::vector<std::string>& soiltypes,
                                     const std::vector<double>& thicknesses,
                                     const double& dz) {
    // ==========================================
    // Biot Gassman equations
    //
    // Parameters:
    // vector[double] KHM: effective bulk modulus over depth [Pa]
    // vector[double] muHM: effective shear modulus over depth [Pa]
    // vector[double] ks: Bulk moduli of grain of each layer [Pa]
    // vector[double] kf: effective compressibility over depth [Pa-1]
    // vector[double] rhob: bulk density over depth [Kg/m3]
    // vector[string] soiltypes: soil types of each layer in the profile
    // vector[double] thicknesses: thickness of each layer [m]
    // double dz: depth discretization [m]
    //
    // Outputs:
    // biotGassmann_Result result: structure with the following fields
    // vector [double] Vp: P-wave velocity over depth [m/s]
    // vector[double] Vs: S-wave velocity over depth [m/s]
    //
    // Programmers:
    // Firstly implemented in Matlab in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // Adaptated for multi-layers by J. Cunha Teixeira in 2024/04
    // ==========================================

    size_t num_points = KHM.size(); // Number of cells in the profile

    std::vector<double> Vp(num_points); // P-wave velocity over depth
    std::vector<double> Vs(num_points); // S-wave velocity over depth

    int start = 0; // Depth start index for fisrt layer

    for (size_t j = 0; j < soiltypes.size(); ++j) {
        std::string_view soilType = soiltypes[j]; // Soil type for current layer
        double thickness = thicknesses[j]; // Thickness of current layer

        selectSoilType_Result soil = selectSoilType_src(soilType); // Soil properties for current layer

        int end = round(thickness / dz) + start; // Depth end index for current layer

        for (int i = start; i < end; ++i) {
            double alpha2, K;
            alpha2 = 1 - KHM[i] / ks[j];
            K = KHM[i] + (alpha2 * alpha2) / (soil.phi / kf[i] + (1 - soil.phi) / ks[j] - KHM[i] / (ks[j] * ks[j]));

            Vp[i] = sqrt((K + (4.0 / 3.0) * muHM[i]) / rhob[i]); // P-wave velocity
            Vs[i] = sqrt(muHM[i] / rhob[i]); // S-wave velocity
        }
        start = end; // Update start index for next layer
    }

    biotGassmann_Result result; // Structure to store the results
    result.Vp = Vp;
    result.Vs = Vs;

    return result;
}



double fish_src(double vp, double vs) {
    // ==========================================
    // Poisson's ration from Vp and Vs
    //
    // Parameters:
    // double vp: P-wave velocity [m/s]
    // double vs: S-wave velocity [m/s]
    //
    // Outputs:
    // double result: Poisson's ratio [-]
    //
    // Programmers:
    // Firstly implemented in Matlab in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // Adaptated for multi-layers by J. Cunha Teixeira in 2024/04
    // ==========================================

    double ratio = vp/vs;
    // return (vp*vp - 2*vs*vs) / (2 * (vp*vp - vs*vs))
    return (0.5 * ratio*ratio - 1) / (ratio*ratio - 1);
}