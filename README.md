![Ginette-2](Ginette-2.png)
==================================================================================
Ginette is a 2-D variably saturated groundwater flow model with integrated 1-D surface flow for the stream. It rigorously simulates water, energy, and solute transport fluxes in porous media using physically-based equations. Ginette was initially developed at METIS (Sorbonne University) to deal with interactions between streams and aquifers as they fluctuate from a connected to a disconnected status. Numerical simulations of experimental laboratory results reproducing such conditions provided the opportunity to test the coupled 1D surface water–2D variably saturated groundwater code (Rivière et al., 2014). Then, Ginette was extended to incorporate coupled heat transfer and water flow in saturated porous media.  The code was compared to experimental data acquired on a complex laboratory system to provide validation on the physical processes and mathematical formulations, in particular for the representation of density change between frozen and liquid water (Rivière et al., 2019). The coupling of fluid flow and heat transfer, accounting for freezing and thawing processes, is implemented in the code for a fully saturated medium and unsaturated zone.

The code is now jointly developed by MINES Paris (PSL University), METIS (Sorbonne University), and Cerege (Marseille University).

 Real-world cryo-hydrogeological paleo-applications, which have been presented in conferences (e.g. Jost, 2011; Jost et al., 2014), were also proposed using Ginette, requiring some additional adaptation to the specific needs of basin-scale calculations. 
The code is applied to estimate the hydraulic and thermal properties of the hyporheic zone (stream-aquifer interface) from the pressure differential and the temperature profle ( Cucchi et al. 2018 and 2021).

Ginette strives to provide the user with complete simulation control. The command files can be created using either a Python or a R script.

## Code organization

The numerical engine is written in Fortran ([src/ginette_V2.f90](src/ginette_V2.f90)), which solves the coupled flow and heat-transport equations described above. It is driven by a Python layer, [src/src_python/](src/src_python/), providing compilation of the Fortran binary, generation of initial and boundary conditions, a generic parameter grid search, comparison to observations, and analytical-solution validation (`Analytical_validation.py`, `Stallman_analysis.py`); this layer is imported by most applications in the `application` directory.

Within `application`, individual studies are implemented in one of two styles:

- **Numbered Python script pipelines**, for reproducible, command-line-driven workflows, typically centered on a single configuration file (e.g., `application/1D_Stream_aquifer_GridSearch`, `application/heat_transport_analytical_validation`, `application/model_dharrma`).
- **Jupyter notebooks**, for interactive exploration, teaching material, or case studies combining simulation and visualization (e.g., `application/RIV2D`, `application/2024_TD_ENS`, `application/2026_TD_ENS`, `application/1D_Diffusive`).

Some applications, such as `application/mini-LOMOS`, combine both Python and R.

Different applications are available in the `application` directory. The most actively maintained ones are described below; each has its own `README.md` with full usage instructions.

### Benchmarking and code verification

**Interfrost benchmark** (`application/Interfrost`)
The Interfrost group proposes a benchmark exercise dealing with subsurface thermo-hydrologic processes involving freezing and thawing, as presented by Painter et al. (2012) and situated within the field of cryohydrogeology (McKenzie et al., 2020). The first phase of the project is limited to Darcy flow in a fully saturated porous medium coupled with heat transfer, advection, and phase change; extensions to the Richards equation and to the air phase are considered for later phases. The benchmark combines test cases inspired by the existing literature (e.g., McKenzie et al., 2007) with new ones, complemented by experimental cases run in a cold room. The project remains open to new or alternative cases of numerical or process-oriented interest to the cold-region community.

**Verification against analytical solutions** (`application/heat_transport_analytical_validation`)
This application verifies the heat transport engine of Ginette against four analytical solutions of reference from the hydrogeological literature, of increasing complexity: Stallman (1965, transient sinusoidal surface temperature signal, homogeneous medium), Bredehoeft and Papadopulos (1965, steady-state temperature profile, homogeneous medium), Shan and Bodvarsson (2004, steady-state profile across a two-layer medium), and Kurylyk et al. (2017, steady-state profile across one to four layers, reproducing their Figure 2). Each case runs Ginette under conditions for which an exact mathematical solution exists and reports the discrepancy between the simulated and the analytical temperature. Across the four cases, Ginette reproduces the reference solutions with errors of a few hundredths to a few tenths of a degree Celsius, including for the most demanding configuration (layered media with strong thermal-conductivity contrasts).

### Unsaturated zone

Several test cases in the unsaturated zone are available (e.g., `application/ZNS-1D`, `application/ZNS_rain`, `application/GINETTE_ZNS`, `application/Warrick`).

### Hyporheic zone characterization from field data (LOMOS-mini)

Two related applications estimate the hydraulic and thermal properties of the hyporheic zone (stream-aquifer interface) from LOMOS-mini sensor data: four temperature probes at 10/20/30/40 cm depth in the streambed together with a pressure-differential sensor (Cucchi et al., 2018, doi:10.1016/j.jhydrol.2017.10.074, hal-01656455).

- `application/mini-LOMOS`: the original inversion pipeline (R and batch-driven), configured through command files (`inversion.COMM`, `inversion_PT100.COMM`, `inversion_parameter.COMM`).
- `application/1D_Stream_aquifer_GridSearch`: a Python re-implementation built around a parameter grid search (permeability, porosity, thermal conductivity, heat capacity), with automated misfit computation (L2, MAE, KGE, PBIAS) and a parameter-identifiability diagnostic based on the thermal Péclet number (Kurylyk et al., 2017), following the approach of Cucchi et al. (2021). It specifically addresses the choice of the upper boundary condition, since the free-water temperature sensor does not measure the same physical quantity as the sediment temperature at the streambed surface.

### Coupled hydro-geophysical modeling

**DHARRMA** (`application/model_dharrma`, Direct HydrogeophysicAl Resistivity and Refraction Modeling Application)
DHARRMA couples the Ginette hydro-thermal model to forward geophysical models for a given infiltration/evaporation scenario and soil facies: electrical resistivity, via Archie's law and the Waxman-Smits model (solved with pyGIMLi), and seismic refraction, via a Hertz-Mindlin rock-physics model (solved with Geopsy). It was developed by N. Radic, A. Rivière, L. Bodet, M. Gautier, A. Gesret, and R. Martin (2025).

### Teaching material

`application/1D_Diffusive` provides a documented notebook illustrating the Ginette energy-transport equation on a simple one-dimensional heat-diffusion case (vertical column, cyclic surface temperature, negligible flow). `application/2024_TD_ENS` and `application/2026_TD_ENS` contain the notebooks (direct model, synthetic cases, MCMC calibration) used for practical sessions at the École Normale Supérieure.


## Authors:
- Riviere, Agnes, agnes.riviere@mines_paristech.fr
- Goncalves, Julio, goncalves@cerege.fr 
- Jost, Anne, anne.jost@upmc.fr


## Contributions by:
- Texier, Jérôme, texier@cerege.fr 


## Userguide
Userguide available in [User_guide.md](User_guide.md)


## References:
- Rivière, A., Gonçalvès, J., Jost, A., & Font, M. (2014). Experimental and numerical assessment of transient stream-aquifer exchange during disconnection. Journal of Hydrology, 517, 574–583. https://doi.org/10.1016/j.jhydrol.2014.05.040
- Rivière, A., Jost, A., Gonçalvès, J. & Font, M. (2018). Pore water pressure evolution below a freezing front under saturated conditions: Large-scale laboratory experiment and numerical investigation. Cold Regions Science and Technology, 158, 76-94. https://doi.org/10.1016/j.coldregions.2018.11.005
- Rivière, A., Gonçalvès, J., Jost, A., (2020). agnes-riviere/ginette: Ginette-2020-09 (Version 2020-09). Zenodo. http://doi.org/10.5281/zenodo.4058821
- Dangeard, M., Rivière, A., Bodet, L., Schneider, S., Guérin, R., Jugnot, D. et Maineult, A. River Corridor Model Constrained by Time‐Lapse Seismic Acquisition. Water Resources Research, American Geophysical Union, 2021, 57 (10), https://doi.org/10.1029/2020WR028911
- Cucchi, K., Flipo, N., Rivière, A. and Rubin, Y. Estimating Hydrothermal Properties and High-Frequency Fluxes From Geophysical Measurements in the Hyporheic Zone. Frontiers in Water, Frontiers, 2021, 3, https://doi.org/10.3389/frwa.2021.700274
- Cucchi, K.,  Rivière, A., Baudin, A., Berrhouma, A., Durand, V., Rejiba, F. et Robin, Y. LOMOS-mini: A coupled system quantifying transient water and heat exchanges in streambeds. Journal of Hydrology, Elsevier, 2018, 561, https:/10.1016/j.jhydrol.2017.10.074


## Message:
"If you use this software, please cite it as below."
Rivière, A., Gonçalvès, J., Jost, A., Ginette,   [![DOI](https://zenodo.org/badge/242535776.svg)](https://zenodo.org/badge/latestdoi/242535776)


authors:
  - Rivière Agnès
    orcid: https://orcid.org/0000-0002-6002-3189
  - Gonçalvès Julio
    orcid: https://orcid.org/0000-0003-0047-4233
  - Jost Anne
    orcid: https://orcid.org/0000-0002-0925-3376
    
    
title: agnes-riviere/ginette: Ginette-2020-09
version: 2020-09
date-released: 2017-12-18

