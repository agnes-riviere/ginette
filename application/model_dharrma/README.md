       /                                                                        \
      /      _____  _    _           _____  _____  __  __                 .      \
     /      |  __ \| |  | |   /\    |  __ \|  __ \|  \/  |   /\          / \      \
    |       | |  | | |__| |  /  \   | |__) | |__) | \  / |  /  \        /   \      |
    |       | |  | |  __  | / /\ \  |  _  /|  _  /| |\/| | / /\ \      |     |     |
    |       | |__| | |  | |/ ____ \ | | \ \| | \ \| |  | |/ ____ \     |     |     |
     \      |_____/|_|  |_/_/    \_\|_|  \_\_|  \_\_|  |_/_/    \_\     \___/     /
      \                                                                          /
       \                              D H A R R M A                             /


# DHARRMA Model (Direct HydrogeophysicAl Resistivity and Refraction Modeling Application)

Code for running the direct transient hydrogeophysical model developed by N. RADIC, A. RIVIERE, L. BODET, S. PASQUET, M. GAUTIER, A. GESRET, R. MARTIN in 2025.
Required inputs: simulation details (days, time steps), facies, soil/thermal/ERT/seismic parameters, infiltration/evaporation scenario, ERT/seismic configuration...

Made up of 6 parts:

    0. Running code section
    1. Model input (infiltration scenario, soil to model, facies, physical parameters, etc.)
    2. Running the hydro and thermal model (see Ginette run)
    3. Running the seismic model (Hertz-Mindlin rock physics model... and forward problem using Géopsy)
    4. Running the electrical model (Archie's/Waxman-Smits law and forward problem using PyGimli)
    5. Data visualization

For initialize the model in a terminal:
    1 Take place in the ginette repository
    2 use the command 'make init_dharrma'

Then tu run the model:
    'make run_dharrma'

## Part 0 : Running code section

You choose witch section of the code will be run. 
- 'lancer_ginette' means that the hydrological simulation will run (True) or not (False). If not other section will use result from the last hydrological simulation
- 'thermique' for thermal simulation. It can run only if lancer_ginette = True
- 'sismique' for seismic model
- 'electrique' for electrical model

## Part 1 : Model input

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
    3. Electrical parameter simulation
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

## Part 2 : Hydro/thermal model

Lunch of the hydrological forward model

## Part 3 : Seismic Model

Lunch of the siesmic forward model

## Part 4 : Electrical Model

Lunch of the Electrical forward model

## Part 5 : Visualisation

All plots are shown


