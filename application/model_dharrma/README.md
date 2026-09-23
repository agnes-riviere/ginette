<img src="logo_dharrma.png" alt="logo_dharrma" width="700">


# DHARRMA Model (Direct HydrogeophysicAl Resistivity and Refraction Modeling Application)

Code for running the direct transient hydrogeophysical model developed by N. RADIC, A. RIVIERE, L. BODET, S. PASQUET, M. GAUTIER, A. GESRET, R. MARTIN in 2025.
Required inputs: simulation details (days, time steps), facies, soil/thermal/ERT/seismic parameters, infiltration/evaporation scenario, ERT/seismic configuration...

The code main_DHARRMA.py is composed of 6 parts:

    0. Running code selection
    1. Model input (infiltration scenario, soil to model, facies, physical parameters, etc.)
    2. Running the hydro and thermal model (see Ginette run)
    3. Running the seismic model (Hertz-Mindlin rock physics model... and forward problem using Géopsy)
    4. Running the electrical model (Archie's/Waxman-Smits law and forward problem using PyGimli)
    5. Data visualization

# For initialize the model in a terminal:
1. Take place in the ginette repository
```bash
cd ginette/
```
2. Use the command :
```bash
make init_dharrma
```

3. Then run the model:
```bash
make run_dharrma
```
If you want to change the simulation parameters or other settings, you can modify main_DHARRMA.py as follows:
### Part 0 : Running code section

You choose witch section of the code will be run. 
- 'lancer_ginette' means that the hydrological simulation will run (True) or not (False). If not other section will use result from the last hydrological simulation
- 'thermique' for thermal simulation. It can run only if lancer_ginette = True
- 'sismique' for seismic model
- 'electrique' for electrical model

### Part 1 : Model input

Definition of all the input you need for the simulation

    1. Ginette Hydro/thermal parameters
        - Mesh
        - Hydrofacies
        - Heterogeneity (work in progress)
        - input file for Ginette. You can create a temperature or a saturation scenario. if Creation_XXX = False simulation use scenario in E_temp_t or E_debis_haut.dat
        - Thermal parameter
    2. Siesmic parameter simulation
        - Start, step and end of simulation (day)
        - Rock physics parameters
        - Geopsy forward model parameters
    3. Electrical parameter simulation (Architecture of the code presented in Solazzi et al. (2021))
        - Start, step and end of simulation (day)
        - Archie parameters
        - Thermal correction
        - Choice of the petrophysical law (archie or Waxmann-smits)
        - pyGIMLy forward model parameters
    4. Visualisation parameter
        - Choice of plot at the end of simulation

    5. (Work in progress) Comparison of simulated WT with a real piezometer
        - need to path to the piezo file (in xlsx)
        - Elevation of the piezo

### Part 2 : Hydro/thermal model

Lunch of the hydrological forward model

### Part 3 : Seismic Model

Lunch of the siesmic forward model
Adaptation of the code developed by Solazzi et al. (2021), extended with a new case: incorporating pressure and saturation computed by a hydrological model.

### Part 4 : Electrical Model

Lunch of the Electrical forward model
Computation of the true resistivity using petrophysical equations, followed by the use of the pyGIMLi library (Rücker et al.) for vertical electrical sounding (VES).

### Part 5 : Visualisation

All plots are shown

## Authors:
- Radic, Nicolas, nicolas.radic@minesparis.psl.eu
- Riviere, Agnes, agnes.riviere@mines_paristech.fr


## References:

- Radic, N., Rivière, A., Bodet, L., Pasquet, S., Martin, R., Gautier, M., Gesret, A. (in prep) : Novel transient hydrogeophysical process-based model tointerpret the geophysical data in the vadose zone. Vadose Zone Journal. This article will be published using this model, and all the figures presented in the article were generated using the figure_article.py script.
- Rücker, C., Günther, T., & Wagner, F. M. (2017). pyGIMLi: An open-source library for modelling and inversion in geophysics. Computers & Geosciences, 109, 106-123. https://doi.org/10.1016/j.cageo.2017.07.011
- Solazzi, S. G., Bodet, L., Holliger, K., & Jougnot, D. (2021). Surface‐wave dispersion in partially saturated soils: The role of capillary forces. Journal of Geophysical Research: Solid Earth, 126(12), e2021JB022074. https://doi.org/10.1029/2021JB022074

## Message:
"If you use this software, please cite it as below."
Radic, N., Rivière, A., Ginette,   [![DOI](URL)


authors:
  - Radic Nicolas
    orcid: https://orcid.org/0009-0009-9275-294X
  - Rivière Agnès
    orcid: https://orcid.org/0000-0002-6002-3189
    
date-released: 2026-09-23
