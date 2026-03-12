#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>

struct vanGen_Result {
    std::vector<double> h, Sw, Swe;
};

struct selectSoilType_Result {
    double wsand, wclay, wsilt, phi, alpha, nvg, theta, Swr;
};

selectSoilType_Result selectSoilType_src(const std::string_view& soiltype);

vanGen_Result vanGen_src(const std::vector<double>& z,
                         const double& WT,
                         const std::vector<std::string>& soiltypes,
                         const std::vector<double>& thicknesses);