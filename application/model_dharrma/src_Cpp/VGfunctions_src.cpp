#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>



struct vanGen_Result{
    std::vector<double> h, Sw, Swe;
};

struct selectSoilType_Result{
    double wsand, wclay, wsilt, phi, alpha, nvg, theta, Swr;
};



selectSoilType_Result selectSoilType_src(const std::string_view& soiltype) {
    // ==========================================
    // Based on:
    // - Carsel, R. F., & Parrish, R. S. (1988). 
    //     Developing joint probability distributions of soil water retention 
    //     characteristics. Water Resource Research, 24(5), 755–769.
    //     https://doi.org/10.1029/wr024i005p00755
    //
    // Parameters:
    // string soiltype: Soil type
    //
    // Outputs:
    // selectSoilType_Result soil: structure with the following fields
    // double wsand: sand ratio [-]
    // double wclay: clay ratio [-]
    // double wsilt: silt ratio [-]
    // double phi: soil porosity [-]
    // double alpha: inverse of entry pressure [1/m]
    // double nvg: Van Genuchten parameter [-]
    // double theta: Van Genuchten parameter [-]
    // double Swr: residual water saturation [-]
    //
    // Programmers:
    // Firstly implemented in Matlab in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // ==========================================

    selectSoilType_Result soil; // Structure to store the results

    if (soiltype == "clay") {
        soil.wsand = 0.149;
        soil.wclay = 0.552;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.38;
        soil.alpha = 0.8;
        soil.nvg = 1.09;
        soil.theta = 0.068;
        soil.Swr = soil.theta / soil.phi;    }
    else if (soiltype == "silt") {
        soil.wsand = 0.058;
        soil.wclay = 0.095;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.46;
        soil.alpha = 1.6;
        soil.nvg = 1.37;
        soil.theta = 0.034;
        soil.Swr = soil.theta / soil.phi;    }
    else if (soiltype == "clayloam") {
        soil.wsand = 0.298;
        soil.wclay = 0.326;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.41; // porosity
        soil.alpha = 1.9;
        soil.nvg = 1.31;
        soil.theta = 0.095;  // residual water content
        soil.Swr = soil.theta / soil.phi;    }
    else if (soiltype == "loam") {
        soil.wsand = 0.4;
        soil.wclay = 0.197;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.43;
        soil.alpha = 3.6;
        soil.nvg = 1.56;
        soil.theta = 0.078;
        soil.Swr = soil.theta / soil.phi;    }
    else if (soiltype == "loamysand") {
        soil.wsand = 0.809;
        soil.wclay = 0.064;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.41;
        soil.alpha = 12.4;
        soil.nvg = 1.28;
        soil.theta = 0.057;
        soil.Swr = soil.theta / soil.phi;    }
    else if (soiltype == "cleansand") {
        soil.wsand = 1;
        soil.wclay = 0;
        soil.wsilt = 0;
        soil.phi = 0.43;
        soil.alpha = 14.5;
        soil.nvg = 2.68;
        soil.theta = 0.045;
        soil.Swr = soil.theta / soil.phi;    }
    else if (soiltype == "sand") {
        soil.wsand = 0.927;
        soil.wclay = 0.029;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.43;
        soil.alpha = 14.5;
        soil.nvg = 2.68;
        soil.theta = 0.045;
        soil.Swr = soil.theta / soil.phi;    }
    else if (soiltype == "sandyclay") {
        soil.wsand = 0.52;
        soil.wclay = 0.43;
        soil.wsilt = 0.05;
        soil.phi = 0.38;
        soil.alpha = 2.7;
        soil.nvg = 1.23;
        soil.theta = 0.1;
        soil.Swr = soil.theta / soil.phi;    }
    else if (soiltype == "sandyclayloam") {
        soil.wsand = 0.543;
        soil.wclay = 0.274;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.43;
        soil.alpha = 1.0;
        soil.nvg = 1.23;
        soil.theta = 0.089;
        soil.Swr = soil.theta / soil.phi;    }
    else if (soiltype == "sandyloam") {
        soil.wsand = 0.634;
        soil.wclay = 0.111;
        soil.wsilt = 1 - soil.wsand - soil.wclay;
        soil.phi = 0.41;
        soil.alpha = 7.5;
        soil.nvg = 1.89;
        soil.theta = 0.065;
        soil.Swr = soil.theta / soil.phi;    }
    else {
        throw std::invalid_argument("Invalid soiltype");
    }

    return soil;
}



vanGen_Result vanGen_src(const std::vector<double>& z,
                         const double& WT,
                         const std::vector<std::string>& soiltypes,
                         const std::vector<double>& thicknesses) {
    // ==========================================
    // Water distribution 
    // Static conditions
    //
    // Parameters:
    // vector[double] z: depths of the profile [m]
    // double WT: depth of the water table [m]
    // vector[string] soiltypes: soil types of each layer in the profile
    // vector[double] thicknesses: thickness of each layer [m]
    //
    // Outputs:
    // vanGen_Result result: structure with the following fields
    // vector [double] h: pressure head over depth [m]
    // vector[double] Sw: total saturation over depth [-]
    // vector[double] Swe: effective wetting phase saturation over depth [-]
    //
    // Programmers:
    // Firstly implemented in Matlab by D. Jugnot in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // Adaptated for multi-layers by J. Cunha Teixeira in 2024/04
    // ==========================================

    size_t num_points = z.size(); // Number of cells in the profile
    double dz = std::abs(z[1] - z[0]); // Depth discretization

    std::vector<double> h(num_points); // Pressure head over depth
    std::vector<double> Swe(num_points); // Effective wetting phase saturation over depth
    std::vector<double> Sw(num_points); // Total saturation over depth

    int start = 0; // Depth start index for fisrt layer

    for (size_t j = 0; j < soiltypes.size(); ++j) {
        std::string_view soilType = soiltypes[j]; // Soil type for current layer
        double thickness = thicknesses[j]; // Thickness of current layer
        
        int end = round(thickness / dz) + start; // Depth end index for current layer

        selectSoilType_Result soil = selectSoilType_src(soilType); // Soil properties for current layer
        double m = 1.0 - 1.0/soil.nvg; // Parameters related to the pore size distribution (Carsel & Parrish, 1988)

        for (int i = start; i < end; ++i) {
            h[i] = z[i] + WT;
            Swe[i] = (h[i] <= 0) ? 1 : std::pow(1 + std::pow(soil.alpha * h[i], soil.nvg), -m);
            Sw[i] = Swe[i] * (1 - soil.Swr) + soil.Swr;
        }
        start = end; // Update start index for next layer
    }

    vanGen_Result result; // Structure to store the results
    result.h = h;
    result.Sw = Sw;
    result.Swe = Swe;

    return result;
}