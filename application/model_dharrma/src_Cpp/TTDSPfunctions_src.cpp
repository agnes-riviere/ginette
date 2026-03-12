#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>



std::string writeVelocityModel_src(const std::vector<double>& thk,
                        const std::vector<double>& vp,
                        const std::vector<double>& vs,
                        const std::vector<double>& rho,
                        const std::string_view& under_layers,
                        const int& n_under_layers) {
    // ==========================================
    // Returns a string of the velocity model in a format expected by GPDC
    //
    // Parameters:
    // vector<double> thk: thickness of each layer [m]
    // vector<double> vp: P-wave velocity of each layer [m/s]
    // vector<double> vs: S-wave velocity of each layer [m/s]
    // vector<double> rho: density of each layer [kg/m^3]
    // string under_layers: Layers to put under the studied soil column on the velocity model in GPDC format
    // int n_under_layers: Number of layers in the under_layers
    //
    // Outputs:
    // string velocity_model: velocity model in GPDC format
    //
    // Programmers:
    // Firstly implemented in Matlab by S. Pasquet in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // ==========================================

    int nl = thk.size(); // Number of layers
    std::stringstream velocity_model; // String to store the velocity model
    velocity_model << nl + n_under_layers << '\n'; // Number of layers in the velocity model

    // Write thickness and velocity in dinver format
    for (int i = 0; i < nl; ++i) {
        if (i == nl - 1 && under_layers.size() == 0) {
            velocity_model << "0 ";
        } else {
            velocity_model << thk[i] << ' ';
        }
        velocity_model << vp[i] << ' ' << vs[i] << ' ' << rho[i] << '\n';
    }

    // Write the under_layers
    if (under_layers.size() > 0) {
        velocity_model << under_layers;
    }

    return velocity_model.str();
}



std::vector<double> firstArrival_src(const std::vector<double>& thk,
                                     const std::vector<double>& vv,
                                     const std::vector<double>& Xdata,
                                     const double& trig) {
    // ==========================================
    // Computes S- or P-wave first arrivals
    //
    // Parameters:
    // vector<double> thk: thickness of each layer [m]
    // vector<double> vv: velocity of each layer [m/s]
    // vector<double> Xdata: distance from the source to the receiver [m]
    // double trig: trigger time [s]
    //
    // Outputs:
    // vector<double> Thod: hodochrone table
    //
    // Programmers:
    // Firstly implemented in Matlab by L. Bodet in Solazzi et al. (2021)
    // Traduced in C++ by J. Cunha Teixeira in 2024/03
    // ==========================================

    size_t rows = thk.size();
    size_t cols = Xdata.size();
    std::vector<double> Tr(rows*cols, 0.0); // arrival time table as a flatten array

    std::vector<double> Thod(cols, 0.0); // hodochrone table
    std::vector<double> Tz(rows, 0.0); // intercepts

    // Calculating intercepts
    for (size_t inl = 0; inl < rows; ++inl) {
        double Tzc = 0.0;
        for (size_t inltz = 0; inltz < inl; ++inltz) {
            Tzc += 2 * thk[inltz] / vv[inltz] * sqrt(1 - vv[inltz] * vv[inltz] / (vv[inl] * vv[inl]));
        }
        Tz[inl] = Tzc;

        for (size_t ix = 0; ix < cols; ++ix) {
            Tr[inl + ix * rows] = Xdata[ix] / vv[inl] + Tz[inl];
        }
    }

    // Finding minimum arrival time for each Xdata point and adding trigger
    for (size_t ix = 0; ix < cols; ++ix) {
        double min_time = Tr[0 + ix * rows];
        for (size_t inl = 1; inl < rows; ++inl) {
            min_time = std::min(min_time, Tr[inl + ix * rows]);
        }
        Thod[ix] = min_time + trig;
    }

    return Thod;
}