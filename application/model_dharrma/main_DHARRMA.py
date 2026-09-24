'''
Modèle DHARRMA (Direct HydrogeophysicAl Resistivity and Refraction Modeling Application)


Code permettant de lancer le modèle direct transitoire hydrogéophysique développé par N. RADIC, A. RIVIERE, L. BODET, S. PASQUET M. GAUTIER, A. GESRET, R. MARTIN en 2025
Input à renseigner : Détail de la simulation (jours, pas), Faciès,parametre sol/thermique/ert/sismique, scénario d'infiltration/evaporation, configuration ERT/sismique...


Constitué de 6 partie :
    0. Paramètre du code
    1. Input du modele (Scénario infiltration, sol à modéliser, faciès, paramètres physique ect. )
    2. Lancement du modèle hydro et thermique (Cf lancement de Ginette)
    3. Lancement du modèle sismique (Modèle rock physic Hertz-mindlin... et du problème direct à l'aide Géopsy)
    4. Lancement du modèle électrique (Loi d'Archie/Waxman-Smits et problème direct à l'aide de Py Gimly)
    5. Visualisation des données

'''

# =====================================================================================
#       /                                                                        \
#      /      _____  _    _           _____  _____  __  __                 .      \
#     /      |  __ \| |  | |   /\    |  __ \|  __ \|  \/  |   /\          / \      \
#    |       | |  | | |__| |  /  \   | |__) | |__) | \  / |  /  \        /   \      |
#    |       | |  | |  __  | / /\ \  |  _  /|  _  /| |\/| | / /\ \      |     |     |
#    |       | |__| | |  | |/ ____ \ | | \ \| | \ \| |  | |/ ____ \     |     |     |
#     \      |_____/|_|  |_/_/    \_\|_|  \_\_|  \_\_|  |_/_/    \_\     \___/     /
#      \                                                                          /
#       \                              D H A R R M A                             /
#
# =====================================================================================


##################### PART 0 : Running code section #########################################

# MUST set matplotlib backend BEFORE importing pyplot
# Use Agg for non-interactive rendering (most robust in all environments)
import matplotlib
matplotlib.use("Agg")

import fonction_DHARRMA as fct
import Direct_model as Dm
import Init_folders as Info
import subprocess
import os
import sys
from matplotlib.cm import copper
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from pathlib import Path
from io import StringIO
from subprocess import run, PIPE, CalledProcessError
import pygimli as pg
from pygimli.physics import VESManager
from lib.VGfunctions import vanGen, selectSoilType
from lib.RPfunctions import hillsAverage, effFluid, hertzMindlin, hertzMindlin_trans, biotGassmann, fish
from lib.TTDSPfunctions import firstArrival, writeVelocityModel, readDispersion

dossier_actuel = Path(__file__).parent


# Selection of code sections to run
lancer_ginette = True
thermique = True
sismique = False
electrique = False
visualisation = True


######################## PART I : MODEL INPUT ##########################################

# General Parameters

nbr_jour =  120 # total number of days to simulate
facies = 'silt' # Facies for the homogeneus simulation (cf Carsel and Parish (1986))


###### GINETTE Hydro Simulation Parameters -----------------------------------------------------------------------------------------------------------------------

# GINETTE Mesh
depth_top = 0 # (m) top of the model
depth_bottom = -4 # (m) bottom of the model
hauteur_WT_initial = -2 # (m) initial water table height
dz = 0.01 # (m) length of the hydrological model cells

# Paramètre hydro
pas_hydro = 900 # (s) Time step use for hydro simulation (also for output files)
homogeneite = True # False = Account for soil heterogeneity (work in progress)
sortie_hauteur_WT = False #True = Create S_wt_depth_t.dat file giving water table height over time


if homogeneite: #  Homogeneus soil parameters
    k1 = 6.94E-14 # Permeability (m2)
    soil = selectSoilType(facies) # Carsel and Parrish's distribution for the selected facies
    phi_soil = soil[3] # Porosity
    Swr_soil = soil[7] # Residual saturation

    #Van genuchten parameter
    alpha_soil = soil[4]
    nvg_soil = soil[5]

else : # Heterogeneous soil parameters (work in progress)
    nbr_couches = 2
    depth_boundary = -3 # (m) depth of the boundary

    #------ Couche 1------
    facies1 = 'sand' #(cf Carsel and Parish (1986))

    k1 = 6.94E-14 # Permeability (m2)
    soil1 = selectSoilType(facies1) # Carsel and Parrish's distribution for the selected facies
    phi_soil1 = soil1[3] # Porosity
    Swr_soil1 = soil1[7] # Residual saturation

    #Van genuchten parameter
    alpha_soil1 = soil1[4]
    nvg_soil1 = soil1[5]

    #------ Couche 2------
    facies2 = 'clay' # Faciès utilisé dans la simulation (cf Carsel and Parish (1986))
    k2 = 1E-14 # Permeability (m2) of the second layer
    soil2 = selectSoilType(facies2) # Carsel and Parrish's distribution for the selected facies
    phi_soil2 = soil2[3] # Porosity
    Swr_soil2 = soil2[7] # Residual saturation

    #Van genuchten parameter
    alpha_soil2 = soil2[4]
    nvg_soil2 = soil2[5]




if thermique:
    ith=1
else:
    ith=0
Dm.setup_ginette_DHARRMA(pas_hydro, ith, nbr_jour, depth_top, depth_bottom, abs(depth_top-depth_bottom), dz, pas_hydro,hauteur_WT_initial)

# Creating input file for GINETTE--------------------------------------------------------------------------------------------------------------------
Creation_temp = False
Creation_infiltation = True
Evapo = False

# Creating the temperature file
# If a temperature file already exists, it can be placed in the E_temp_t.dat file
if Creation_temp :

    bottom_temp = 10 #(°C) stable temp at the bottom of the model
    deg_per_day = 0.2 #0#(°C) Linear temperature increase over time
    daily_fluctuation = 5
    weekly_fluctuation = 3
    temp_offset = 0

    fct.creation_temp(nbr_jour,pas_hydro,bottom_temp,deg_per_day,daily_fluctuation,weekly_fluctuation,temp_offset)

# Creating the infiltration file
# If an infiltration file already exists, it can be placed in the E_debit_haut_t.dat file
if Creation_infiltation :
    min_infiltration = 0
    max_infiltration = 4E-08 # m/s
    duree_transition = 1#  Days
    duree_max_infiltration= 1 # Days
    debut_jour_infiltration = [7,12,17,22,30,87,92,97,102] #[30,102] #Starting days of infiltration events


    if Evapo :
        max_evapo = 7E-8
        duree_transition_evapo = 0.125 # Days
        duree_max_evapo = 0.5 # Days
        debut_jour_evapo = [i+0.5 for i in range(1,10)]
        fct.creation_infiltration_evapo(nbr_jour,pas_hydro,min_infiltration,max_infiltration,duree_transition,debut_jour_infiltration,duree_max_infiltration,0,max_evapo,duree_transition_evapo,debut_jour_evapo,duree_max_evapo)
    else :
        fct.creation_infiltration(nbr_jour,pas_hydro,min_infiltration,max_infiltration,duree_transition,debut_jour_infiltration,duree_max_infiltration)

###### Paramètre de simulation Thermique GINETTE--------------------------------------------------------------------------------------------------------------------

lambda1 = 2.4
C1 = 897
rho1 = 2400

if homogeneite == False:
    lambda2 = 1.5
    C2 = 800
    rho2 = 2000


if homogeneite:
    Dm.generate_zone_parameters_DHARRMA(thermique,depth_bottom,depth_top, dz, k1, phi_soil, alpha_soil, nvg_soil, Swr_soil, lambda1, C1, rho1)
else :
    Dm.generate_zone_parameters_hetero_DHARRMA(thermique,depth_bottom,depth_top, dz,nbr_couches, depth_boundary, k1, phi_soil1, alpha_soil1, nvg_soil1, Swr_soil1, lambda1, C1, rho1,
                                        k2, phi_soil2, alpha_soil2, nvg_soil2, Swr_soil2, lambda2, C2, rho2)


###### Seismic simulation (with Geopsy) --------------------------------------------------------------------------------------------------------------
######### Parametre simulation #############
debut_sim_sis = 1 # (day) start of the simulation
pas_sim_sis = 1 # (day) step
fin_sim_sis = 120#(day inclued) end

dz_sis = 0.01 # (m) vertical discretization

first_arrival_calcul = False

######### ROCK PHYSICS PARAMETERS ##########
# General physical constants
rhow = 1000.0 # Water density [Kg/m3]
rhoa = 1.0 # Air density [Kg/m3]
kw = 2.3e9 # Water bulk modulus [Pa]
ka = 1.01e5 # Air bulk modulus [Pa]
g = 9.82 # Gravity acceleration [m/s2]

# Grains/agregate mechanical properties
mu_clay = 6.8 # Shear moduli [GPa]
mu_silt = 45.0
mu_sand = 45.0
k_clay = 25.0 # Bulk moduli [GPa]
k_silt = 37.0 
k_sand = 37.0 
rho_clay = 2580.0 # Density [kg/m3]
rho_silt = 2600.0
rho_sand = 2600.0

if homogeneite:
    soiltypes = [facies]
else:
    soiltypes = [facies1,facies2]
    # soiltypes2 = [facies2]

# Grains/agregate parameters per layer
if homogeneite:
    Ns = [9] # Coordination Number (number of contact per grain) | default = 8
    fracs = [0.3] # Fraction of non-slipping grains (helps making the soil less stiff) | default = 0.3
else:
    Ns = [9,9] # Coordination Number (number of contact per grain) | default = 8
    fracs = [0.3,0.3] # Fraction of non-slipping grains (helps making the soil less stiff) | default = 0.3
# Four possible RP models:
# kk = 1 # Constant Pe (see the approach of Zyserman et al., 2017)
# kk = 2 # Pe without suction
# kk = 3 # Pe with suction (cf. Solazzi et al. 2021)
kk = 4 # External pression data (From ginette output)

######### Parametre modèle direct geopsy ##############
# In GPDC format : [thickness Vp Vs rho]
under_layers = [] # Empty list if no under layers
# under_layers = [
#                 [10, 4000, 2000, 2500],
#                 [0, 8000, 4000, 2500],
#                 ]
under_layers = np.array(under_layers)
n_under_layers = len(under_layers) # Number of under layers

x0 = 0.125 # first geophone position [m]
Nx = 96 # number of geophones [m]
dx = 0.25 # geophone interval [m]
xs = np.arange(x0, Nx * dx + 1, dx)
trig  = 0 # data pretrig (if needed)


# Frequency domain and sampling setup to compute dispersion
nf = 500 # number of frequency samples [#]
df = 1 # frequency sample interval [Hz]
min_f = 15 # minimum frequency [Hz]
max_f = min_f + (nf - 1) * df


n_modes = 1 # Number of modes to compute
s = 'frequency' # Over frequencies mode
wave = 'R' # Rayleigh (PSV) fundamental mode


###### Electrical simulation(with PyGImly) -----------------------------------------------------------------------------------------------------------
######### Parametre simulation #############
debut_sim_elec = 1 # (day) start of the simulation
pas_sim_elec = 1 # (day) step
fin_sim_elec = 120#(day inclued) end

elec_static = True # Simulation of an hydrostatic model

### Petrophysical parameters/law

#Parametre Loi d'Archie
a_archie = 1.196 # Facteur de tortuosité ]0.5;1.5] 
m_archie = 1.929 # Exposant de concentration
n_archie = 2.338 # Exposant de saturation
# Parametre a,m et n calculé avec les techniques conventionelles.

if homogeneite ==  False:
     a_archie2 = 1.147 # Facteur de tortuosité ]0.5;1.5] #limon = 1.196, sable = 1.147
     m_archie2 = 2.135 # Exposant de concentration # limon = 1.929, sable = 2.135
     n_archie2 = 0.858 # Exposant de saturation # limon = 2.338, sable = 0.858

beta_s = 5.2E-9

#Thermal correction 
rho_water_25 = 75 # Resistivité elec du fluide ici de l'eau à 25°C (en ohm.m)
a_T = 0.02 # Compensateur de Température Hayashi 2004 et valeur de Matthes 1982

Waxman_smits = True
if Waxman_smits:
    if homogeneite:
        wsand = soil[0]
        wclay = soil[1]
        wsilt = soil[2]
        CEC = wclay*19270 # illite Woodruff and Revil 2011
        Q_v = rho1*((1-phi_soil)/phi_soil)*CEC

        print("Q_v =",Q_v)
    else :
        wsand1 = soil1[0]
        wclay1 = soil1[1]
        wsilt1 = soil1[2]
        CEC1 = wclay1*19270 # illite Woodruff and Revil 2011
        Q_v1 = rho1*((1-phi_soil1)/phi_soil1)*CEC1

        wsand2 = soil2[0]
        wclay2 = soil2[1]
        wsilt2 = soil2[2]
        CEC2 = wclay2*19270 # illite Woodruff and Revil 2011
        Q_v2 = rho2*((1-phi_soil2)/phi_soil2)*CEC2

        
        print("Q_v1 =",Q_v1)
        print("Q_v2 =",Q_v2)


### FWD model : apparent resistivity

ab2_scale = 'log' # Recommanded
# ab2_scale = 'linear' 
ab2_min = 1.5 # m
ab2_max = 100 # m
ab2_nbr_pt = 100 # Nombre d'écartement des électrodes AB

mn2 = 1.0 # m, écartement des électrodes MN

###### Paramètre Visualisation -------------------------------------------------------------------------------------------------------------------------------------

visualisation_temp = True   # Visualisation of the temperature scenario
visualisation_pluie = True # Visualisation of the infiltration scenario

#### Visualisation 2D ####
debut_representation = 1 # (day) Start of the representation
fin_representation = nbr_jour # (day) End of the representation
pas_representation = 1 # (day) Step of the representation

lim_depth = 2.1 # (m) Depth up to which the soil is represented
#ginette output 2D
visualisation_output_ginette = True # Visualisation of hydrological parameters in 2D

#FWD model
visualisation_propriete_geophy_2D = False   # Visualisation of geophysical properties in 2D
visualisation_observable_geophy_2D = False # Visualisation of geophysical observables in 2D

#Daily profile representation
jour_profil = [1,30,40,60] # Day of the profile representation
representation = 1 # All representations on the same graph
# representation = 2 # One representation per graph
visualisation_propriete_hydro_profil = True # Profile of hydrological properties (saturation, pressure, temperature)
visualisation_propriete_geophy_profil = True # Profile of geophysical properties (Vp, Vs, rho)
visualisation_observable_geophy_profil = True # Profile of geophysical observables (dispersion curves, apparent resistivity)
visualisation_wt = False

#Comparaison WT simulé avec un vrai piezomètre (work in progress)

comparaison_wt_piezo = False
date_debut_simulation = "01/01/2025 00:00"  # Format: JJ/MM/AAAA HH:MM
cote_ngf_piezo = 80.02 # Cote en mNGF du piezo que l'on va comparer
piezo_data_dir = "/home/ariviere/Documents/Bassin_Orgeval/Hydro_data/processed_data/AvAv2"  # Chemin du repertoire piezometre
piezo_name = "pzps16"  # Nom du piezometre (sans extension)

path_piezo = None
if comparaison_wt_piezo:
    path_piezo = fct.select_piezometer_file(piezo_data_dir, piezo_name)

DEBUG = False

#################### PARTIE II : MODÈLE HYDRO/THERMIQUE ######################################
depth = depth_top-depth_bottom
zs = -np.arange(dz, depth + dz, dz) # Depth positions (negative downward) [m]
thks = np.diff(np.abs(zs)) # thickness vector [m]


# Compilation Ginette
Info.compile_ginette_DHARRMA(DEBUG)

if lancer_ginette: 
    # Lancement Ginette
    # Changez de répertoire pour build
    os.chdir('input_ginette')
    # fichier = "E_zone.dat"

    # if os.path.isfile(fichier):
    #     print(f"Le fichier '{fichier}' existe et est accessible.")
    # else:
    #     print(f"Le fichier '{fichier}' n'existe pas ou n'est pas accessible.")
    subprocess.call(["./ginette"]) # Utilisation Linux
    # subprocess.call(["ginette.exe"]) # Utilisation Windows
    # Revenez au répertoire précédent
    os.chdir('..')
    print('Run model hydro et thermique Ok')


    if sortie_hauteur_WT:
        # Récupération de la hauteur de la nappe
        fct.creation_S_wt_depth(zs)

####################### PARTIE III : MODÈLE SISMIQUE #########################################

if sismique:

    # Modèle de physique des roches
    pressure = pd.read_csv("input_ginette/S_pressure_profil_t.dat", header=None, sep=r'\s+', names=['dt', 'Z', 'Pr','h'])
    saturation = pd.read_csv("input_ginette/S_saturation_profil_t.dat", header=None, sep=r'\s+', names=['dt', 'Z', 'Sw'])

    interval_sim_sis = np.arange(debut_sim_sis*86400,fin_sim_sis*86400+pas_sim_sis*86400,pas_sim_sis*86400)
    color_map = copper(np.linspace(0, 1, len(interval_sim_sis)))  # Colorscale for plots if several times are tested

    # print(saturation)
    result = saturation[saturation["dt"].isin(interval_sim_sis)]
    dict_result = result.groupby("dt")["Sw"].apply(list).to_dict()

    it_tot_simique = int((fin_sim_sis-debut_sim_sis)/pas_sim_sis)
    for i, temps in enumerate(dict_result):


        saturation_profil = dict_result[temps]
        pression_profil = pressure.loc[pressure['dt'] == temps, 'Pr'].tolist()
        pression_profil_hPa = [i/1000 for i in pression_profil]

        #Trouver la hauteur de la WT
        z_sat = None
        j=0
        while (j < len(zs) and z_sat == None) :
            if saturation_profil[j] == 1.0 :
                z_sat = zs[j]
            else :
                j = j+1
        #Créer le fichier hs qui correspond à la hauteur au dessus de la WT
        hs = list()
        Swes = list()
        for profondeur in zs : # Création de hs
            h_val = profondeur-z_sat
            
            
            hs.append(round(float(h_val),3))
        
        if homogeneite:
            for prof in range (len(hs)):
                Swes_val = (saturation_profil[prof] - Swr_soil) / (1 - Swr_soil)

                Swes.append(Swes_val)
            thicknesses = [depth]
        else:
            for prof in range (len(hs)):
                if prof <= int(abs(depth_boundary/dz)) :
                    Swes_val = (saturation_profil[prof] - Swr_soil1) / (1 - Swr_soil1)
                else :
                    Swes_val = (saturation_profil[prof] - Swr_soil2) / (1 - Swr_soil2)

                Swes.append(Swes_val)
            thicknesses = [abs(depth_boundary),abs(depth_top-depth_bottom)-abs(depth_boundary)]
        if len(set(map(len, (soiltypes, thicknesses, Ns, fracs)))) != 1:
            raise ValueError(f"Arrays are not the same size : {soiltypes = }, {thicknesses = }, {Ns = }, {fracs = }")


        # Effective Grain Properties (constant with depth)
        mus, ks, rhos, nus = hillsAverage(mu_clay, mu_silt, mu_sand, rho_clay,
                                                rho_silt, rho_sand, k_clay, k_silt,
                                                k_sand, soiltypes)
        

        # Effective Fluid Properties
        kfs, rhofs, rhobs = effFluid(saturation_profil, kw, ka, rhow,
                                        rhoa, rhos, soiltypes, thicknesses , dz) # Utilisation de la saturation total


        # Hertz Mindlin Frame Properties
        KHMs, muHMs = hertzMindlin_trans(Swes, zs, hs, rhobs, pression_profil,
                                    g, rhoa, rhow, Ns,
                                    mus, nus, fracs, kk,
                                    soiltypes, thicknesses) # Utilisation de la saturation effective

        # Saturated Properties
        VPs, VSs = biotGassmann(KHMs, muHMs, ks, kfs,
                                    rhobs, soiltypes, thicknesses, dz)
            

    # SEISMIC FWD MODELING -----------------------------------------------------------------------------------------------------------------------

        # First arrival time computations
        thks_tmp = np.copy(thks) # Thicknesses of the layers
        VPs_tmp = np.copy(VPs) # P-wave velocities of the layers
        VSs_tmp = np.copy(VSs) # S-wave velocities of the layers
        for layer in under_layers:
            thickness = layer[0]
            if thickness == 0:
                thickness = 2*dz
            vp = layer[1]
            vs = layer[2]
            thks_tmp = np.concatenate((thks_tmp, [dz]*int(thickness/dz))) # Thicknesses of the layers
            VPs_tmp = np.concatenate((VPs_tmp, [vp]*int(thickness/dz))) # P-wave velocities of the layers
            VSs_tmp = np.concatenate((VSs_tmp, [vs]*int(thickness/dz))) # S-wave velocities of the layers            
        if first_arrival_calcul:
            ThodPs = firstArrival(thks_tmp, VPs_tmp, xs, trig) # P-wave first arrival times
            ThodSs = firstArrival(thks_tmp, VSs_tmp, xs, trig) # S-wave first arrival times

            # print(ThodPs)
        

        # Velocity model in string format for GPDC
        under_layers_str = '\n'.join([' '.join(map(str, layer)) for layer in under_layers]) + '\n'
        velocity_model_string = writeVelocityModel(thks, VPs, VSs, rhobs, under_layers_str, n_under_layers)

        # Dispersion curves computing with GPDC
        velocity_model_RAMfile = StringIO(velocity_model_string) # Keep velocity model string in the RAM in a file format alike to trick GPDC which expects a file
        gpdc_command = [f"/usr/local/Geopsy.org/bin/gpdc -{wave} {n_modes} -n {nf} -min {min_f} -max {max_f} -s {s}"]

        try:
            process = run(gpdc_command, input=velocity_model_RAMfile.getvalue(), text=True, shell=True, stdout=PIPE, stderr=PIPE, check=True) # Raw output string from GPDC
        except CalledProcessError as e:
            print(f"\nERROR during GPDC computation. Returned:\n{e.stdout}")
            print("Used parameters:")
            print(f'{soiltypes = }')
            print(f'{thicknesses = }')
            print(f'{Ns = }')
            print(f'{fracs = }')
            print(f'{z_sat = }')
            print(f'{dz = }\n')
            print('INFO : Try to reduce dz\n')
            raise

        gpdc_output_string = process.stdout # Raw output string from GPDC
        dispersion_data, n_modes = readDispersion(gpdc_output_string) # Reads GPDC output and converts dispersion data to a list of numpy arrays for each mode
        print(f'Simulation sismique : {i+1}/{it_tot_simique+1}')
        ### SAUVEGARDE DES DONNEES #####
        if i == 0 :
            if homogeneite:
                path_vp_vs = f'sismique/{soiltypes[0]}/output_SL_kk{kk}_Vp_Vs.dat'
                path_PS_v = f'sismique/{soiltypes[0]}/output_SL_kk{kk}_PS_v_Phase.dat'
            
            else:
                path_vp_vs = f'sismique/{soiltypes[0]}_{soiltypes[1]}/output_SL_kk{kk}_Vp_Vs.dat'
                path_PS_v = f'sismique/{soiltypes[0]}_{soiltypes[1]}/output_SL_kk{kk}_PS_v_Phase.dat'

            if not os.path.exists(path_vp_vs):
                os.makedirs(os.path.dirname(path_vp_vs), exist_ok=True)
                with open(path_vp_vs, "w") as f:
                    pass
            if not os.path.exists(path_PS_v):
                os.makedirs(os.path.dirname(path_PS_v), exist_ok=True)
                with open(path_PS_v, "w") as f:
                    pass
            
            if first_arrival_calcul:
                path_first_arrival = f'sismique/{soiltypes[0]}/output_SL_kk{kk}_1AP_1AV_Phase.dat'

                if not os.path.exists(path_first_arrival):
                    os.makedirs(os.path.dirname(path_first_arrival), exist_ok=True)
                    with open(path_first_arrival, "w") as f:
                        pass
                
                with open(path_first_arrival,'w') as f:
                    for xs_data,first_arr_p,first_arr_v in zip(xs,ThodPs,ThodSs):
                        f.write(f"{temps} {xs_data} {first_arr_p}  {first_arr_v}\n")

            with open(path_vp_vs,'w') as f:
                for z_value,Vp_value,Vs_value,rhob in zip(zs,VPs,VSs,rhobs):
                    f.write(f"{temps} {z_value} {Vp_value}  {Vs_value}  {rhob}\n")
            
            with open(path_PS_v,'w') as f:
                for f_value,Vp_value in zip(dispersion_data[0][:,0], dispersion_data[0][:,1]):
                    f.write(f"{temps} {f_value} {Vp_value}\n")
            

        else :
             with open(path_vp_vs,'a') as f:
                for z_value,Vp_value,Vs_value,rhob in zip(zs,VPs,VSs,rhobs):
                    f.write(f"{temps} {z_value} {Vp_value}  {Vs_value}  {rhob}\n")        

             with open(path_PS_v,'a') as f:
                for f_value,PS_v_value in zip(dispersion_data[0][:,0], dispersion_data[0][:,1]):
                    f.write(f"{temps} {f_value} {PS_v_value}\n")

             if first_arrival_calcul:
                 with open(path_first_arrival,'a') as f:
                    for xs_data,first_arr_p,first_arr_v in zip(xs,ThodPs,ThodSs):
                        f.write(f"{temps} {xs_data} {first_arr_p}  {first_arr_v}\n")
                

    print('Sismique OK')

####################### PARTIE IV : MODÈLE ELECTRIQUE ########################################

if electrique:
    depth = depth_top-depth_bottom
    zs = -np.arange(dz, depth + dz, dz) # Depth positions (negative downward) [m]
    thks = np.diff(np.abs(zs)) # thickness vector [m]

    saturation = pd.read_csv("input_ginette/S_saturation_profil_t.dat", header=None, sep=r'\s+', names=['dt', 'Z', 'Sw'])
    temp = pd.read_csv("input_ginette/S_temperature_t.dat", header=None, sep='\s+', names=['dt', 'Z', 'temp'])

    interval_sim_elec = np.arange(debut_sim_elec*86400,fin_sim_elec*86400+pas_sim_elec*86400,pas_sim_elec*86400)

    result = saturation[saturation["dt"].isin(interval_sim_elec)]

    dict_result = result.groupby("dt")["Sw"].apply(list).to_dict()
    it_tot_elec = int((fin_sim_elec-debut_sim_elec)/pas_sim_elec)
    for i, temps in enumerate(dict_result):
        print(f'Simulation Electrique : {i+1}/{it_tot_elec+1}')
        temperature_profil = temp[temp["dt"].isin([temps])]
        temperature = np.array(temperature_profil['temp'])
        # print(np.array(temperature))
        saturation_profil = dict_result[temps]
        sat_profil_array = np.array(saturation_profil)

        rho_water_T = (rho_water_25)/(1+a_T*(temperature - 25)) # Matthes 1982

        if Waxman_smits :
            # kappa = a_archie*(pow(sat_profil_array,n_archie)/pow(phi_soil,-m_archie))*(1/rho_water_T+kappa_s/sat_profil_array)
            # kappa = (pow(sat_profil_array,n_archie)/(a_archie*pow(phi_soil,-m_archie)))*(1/rho_water_T+kappa_s/sat_profil_array)
            # rho_vrai = 1/kappa

            if homogeneite :
                kappa_s = beta_s*Q_v/sat_profil_array
                kappa = (pow(sat_profil_array,n_archie)/(a_archie*pow(phi_soil,-m_archie)))*((1/rho_water_T)+kappa_s)
            
            else :
                 # Creation array Q_v
                kappa = np.zeros(abs(int((depth_bottom-depth_top)/dz)))

                for maille in range (abs(int((depth_bottom-depth_top)/dz))):
                    if maille <= abs(int ((depth_boundary-depth_top)/dz)) : # Test pour savoir dans quelle couche on est
                        
                        #Couche 1
                        kappa[maille] = (pow(sat_profil_array[maille],n_archie)/
                                         (a_archie*pow(phi_soil1,-m_archie)))*((1/rho_water_T[maille])+(beta_s*Q_v1/sat_profil_array[maille]))

                    else :

                        #Couche 2
                        kappa[maille] = (pow(sat_profil_array[maille],n_archie2)/
                                         (a_archie*pow(phi_soil2,-m_archie2)))*((1/rho_water_T[maille])+(beta_s*Q_v2/sat_profil_array[maille]))

            rho_vrai = 1/kappa
            
            # rho_vrai = a_archie*rho_water_T*pow(phi_soil,-m_archie)*pow(sat_profil_array,-n_archie)*((1+rho_water_T*B_WS*Qv_WS)/(1+(rho_water_T*B_WS*Qv_WS/sat_profil_array)))
            # kappa = 5.2E-9


        else : # Loi d'Archie 
            if homogeneite :
                rho_vrai = a_archie*rho_water_T*pow(phi_soil,-m_archie)*pow(sat_profil_array,-n_archie)
            else :
                raise NotImplementedError("Loi d'Archie non implémentée pour le cas hétérogène")

    #### Problème direct (Utilisation de PyGimly)
    # Espacement AB/2
        if ab2_scale == 'log':
            ab2 = np.logspace(np.log10(ab2_min), np.log10(ab2_max), ab2_nbr_pt)
        elif ab2_scale == 'linear':
            ab2 = np.linspace(ab2_min,ab2_max,ab2_nbr_pt)

        # Modèle de rho calculé
        synthModel = pg.cat([dz for i in range(len(zs)-1)], [rho_vrai[j] for j in range(len(rho_vrai))]) # Epaisseur de couche, résistivité correspondante (+ socle)
        # print(len([dz for i in range(len(zs)-1)]))
        ves = VESManager()
        ra = ves.simulate(synthModel, ab2=ab2, mn2=mn2)#, noiselevel = 0.01)
        u = ves.simulate(synthModel, ab2=ab2, mn2=mn2,returnFields = True)# A B M N
        # print(len(u))
    #### SAUVEGARDE DES DONNEES #####
        if i == 0 :
            if homogeneite :
                path_rho_vrai = f'elec/{soiltypes[0]}/rho_vrai.dat'
                path_rho_app_AB2 = f'elec/{soiltypes[0]}/rho_app_AB2.dat'
            else :
                path_rho_vrai = f'elec/{facies1}_{facies2}/rho_vrai.dat'
                path_rho_app_AB2 = f'elec/{facies1}_{facies2}/rho_app_AB2.dat'
            if not os.path.exists(path_rho_vrai):
                os.makedirs(os.path.dirname(path_rho_vrai), exist_ok=True)
                with open(path_rho_vrai, "w") as f:
                    pass
            if not os.path.exists(path_rho_app_AB2):
                os.makedirs(os.path.dirname(path_rho_app_AB2), exist_ok=True)
                with open(path_rho_app_AB2, "w") as f:
                    pass

            with open(path_rho_vrai,'w') as f:
                for z_value,rho_vrai_value in zip(zs,rho_vrai):
                    f.write(f"{temps} {z_value} {rho_vrai_value}\n")
            
            with open(path_rho_app_AB2,'w') as f:
                for ab2_value,ra_value in zip(ab2,ra):
                    f.write(f"{temps} {ab2_value} {ra_value}\n")
        else :
            with open(path_rho_vrai,'a') as f:
                for z_value,rho_vrai_value in zip(zs,rho_vrai):
                    f.write(f"{temps} {z_value} {rho_vrai_value}\n")  

            with open(path_rho_app_AB2,'a') as f:
                for ab2_value,ra_value in zip(ab2,ra):
                    f.write(f"{temps} {ab2_value} {ra_value}\n")

    print('Electrique OK')



if elec_static:
    depth = depth_top-depth_bottom
    thicknesses = [depth]
    zs = -np.arange(dz, depth + dz, dz) # Depth positions (negative downward) [m]
    thks = np.diff(np.abs(zs)) # thickness vector [m]
    # path_saturation =  '/home/nradic/Documents/ginette/application/model_dharrma/save_data/data_article_1/scenario_5_pluies/input_ginette/S_saturation_profil_t.dat'
    path_saturation = 'input_ginette/S_saturation_profil_t.dat'
    # path_temp = '/home/nradic/Documents/ginette/application/model_dharrma/save_data/data_article_1/scenario_5_pluies/input_ginette/S_temperature_t.dat'
    path_temp = 'input_ginette/S_temperature_t.dat'
    saturation = pd.read_csv(path_saturation, header=None, sep='\s+', names=['dt', 'Z', 'Sw'])
    temperature = pd.read_csv(path_temp, header=None, sep='\s+', names=['dt', 'Z', 'temp'])

    for i,jours in enumerate(jour_profil):
        temps_ginette = jours*86400

        saturation_jours = saturation[saturation["dt"]==temps_ginette]
        sat = np.array(saturation_jours['Sw'])
        temp_jours = temperature[temperature['dt']==temps_ginette]
        temp = np.array(temp_jours['temp'])

        # Determination de la hauteur de la WT
        z_sat = None
        j=0
        while (j < len(zs) and z_sat == None) :
            if sat[j] == 1.0 :
                z_sat = zs[j]
            else :
                j = j+1

        # Calcule du profile de saturation en permanent à l'aide de la hauteur de nappe
        hs, Sws, Swes = vanGen(zs, -z_sat, soiltypes, thicknesses)
        Sws_array = np.array(Sws)

        # path_sat_static = '/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/saturation_profil.dat'
        path_sat_static = 'hydrostatic_model/hydro/saturation_profil.dat'
        if i == 0 : # Ecriture dans un fichier du profile de sat permanent
            if not os.path.exists(path_sat_static):
                os.makedirs(os.path.dirname(path_sat_static), exist_ok=True)
            with open(path_sat_static,'w') as f:
                for z_value,Sws_value in zip(zs,Sws):
                    f.write(f"{temps_ginette} {z_value} {Sws_value}\n")
        else :
             with open(path_sat_static,'a') as f:
                for z_value,Sws_value in zip(zs,Sws):
                    f.write(f"{temps_ginette} {z_value} {Sws_value}\n")  

        
        rho_water_T = (rho_water_25)/(1+a_T*(temp - 25))

        # Loi d'Archie ou Waxman-Smits
        if Waxman_smits :

            kappa_s = beta_s*Q_v/Sws_array
            kappa = (pow(Sws_array,n_archie)/(a_archie*pow(phi_soil,-m_archie)))*((1/rho_water_T)+kappa_s)
            # kappa = (pow(Sws_array,n_archie)/(a_archie*pow(phi_soil,-m_archie)))*((1/25)+0.5)
            rho_vrai = 1/kappa
        
        else : # Loi d'Archie 
            rho_vrai = a_archie*rho_water_T*pow(phi_soil,-m_archie)*pow(sat_profil_array,-n_archie)


        # if i == 0 : # Ecriture rho vrai permanent
        #     if not os.path.exists('/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_vrai_static_temperature.dat'):
        #         os.makedirs(os.path.dirname('/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_vrai_static_temperature.dat'), exist_ok=True)
        #     with open(f'/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_vrai_static_temperature.dat','w') as f:
        #         for z_value,rho_vrai_value in zip(zs,rho_vrai):
        #             f.write(f"{temps_ginette} {z_value} {rho_vrai_value}\n")
        # else :
        #      with open(f'/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_vrai_static_temperature.dat','a') as f:
        #         for z_value,rho_vrai_value in zip(zs,rho_vrai):
        #             f.write(f"{temps_ginette} {z_value} {rho_vrai_value}\n")  
        if i == 0 : # Ecriture rho vrai permanent
            if not os.path.exists('hydrostatic_model/elec/rho_vrai_static_temperature.dat'):
                os.makedirs(os.path.dirname('hydrostatic_model/elec/rho_vrai_static_temperature.dat'), exist_ok=True)
            with open(f'hydrostatic_model/elec/rho_vrai_static_temperature.dat','w') as f:
                for z_value,rho_vrai_value in zip(zs,rho_vrai):
                    f.write(f"{temps_ginette} {z_value} {rho_vrai_value}\n")
        else :
             with open(f'hydrostatic_model/elec/rho_vrai_static_temperature.dat','a') as f:
                for z_value,rho_vrai_value in zip(zs,rho_vrai):
                    f.write(f"{temps_ginette} {z_value} {rho_vrai_value}\n")  
    
        # FWD MODEL -----------------------------------------------------------------------------------------------------------------------------------

        # Espacement AB/2 
        ab2 = np.logspace(np.log10(1.5), np.log10(100), 100)
        mn2 = 1.0
        # Modèle de rho calculé
        synthModel = pg.cat([dz for i in range(len(zs)-1)], [rho_vrai[j] for j in range(len(rho_vrai))])
        ves = VESManager()
        ra = ves.simulate(synthModel, ab2=ab2, mn2=mn2)

        # if i == 0 : # Ecriture rho mesuré permanent
        #     if not os.path.exists('/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_app_AB2_static_temperature.dat'):
        #         os.makedirs(os.path.dirname('/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_app_AB2_static_temperature.dat'), exist_ok=True)
        #     with open(f'/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_app_AB2_static_temperature.dat','w') as f:
        #         for ab2_value,ra_value in zip(ab2,ra):
        #             f.write(f"{temps_ginette} {ab2_value} {ra_value}\n")
        # else :
        #      with open(f'/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_app_AB2_static_temperature.dat','a') as f:
        #         for ab2_value,ra_value in zip(ab2,ra):
        #             f.write(f"{temps_ginette} {ab2_value} {ra_value}\n")
        if i == 0 : # Ecriture rho mesuré permanent
            if not os.path.exists('hydrostatic_model/elec/rho_app_AB2_static_temperature.dat'):
                os.makedirs(os.path.dirname('hydrostatic_model/elec/rho_app_AB2_static_temperature.dat'), exist_ok=True)
            with open(f'hydrostatic_model/elec/rho_app_AB2_static_temperature.dat','w') as f:
                for ab2_value,ra_value in zip(ab2,ra):
                    f.write(f"{temps_ginette} {ab2_value} {ra_value}\n")
        else :
             with open(f'hydrostatic_model/elec/rho_app_AB2_static_temperature.dat','a') as f:
                for ab2_value,ra_value in zip(ab2,ra):
                    f.write(f"{temps_ginette} {ab2_value} {ra_value}\n")

################### PARTIE V : VISUALISATION DES DONNÉES #####################################

if visualisation :
    # Create output directory for figures
    from datetime import datetime
    output_dir = dossier_actuel / f"output_figures_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    output_dir.mkdir(exist_ok=True)
    print(f"Figures will be saved to: {output_dir}")
    
    if homogeneite:
        geophy_folder = facies
    else:
        geophy_folder = f"{facies1}_{facies2}"
    required_vp_vs_file = dossier_actuel / "sismique" / geophy_folder / f"output_SL_kk{kk}_Vp_Vs.dat"

    if visualisation_temp :
        fct.plot_temp_scenario(f'{dossier_actuel}/input_ginette/E_temp_t.dat')
    if visualisation_pluie:
        fct.plot_graphe_pluie(f'{dossier_actuel}/input_ginette/E_debit_haut_t.dat',nbr_jour, pas_hydro,path_pluie_B = None, 
                              liste_mesure = jour_profil, hauteur_WT_A = None, hauteur_WT_B = None, barre_vertical = [],path_fig =str(output_dir / "plot_pluie.png"), max_cumul=1E-3)
    if visualisation_output_ginette:
        fct.three_plot_output_ginette_dharrma(dossier_actuel,debut_representation,fin_representation,pas_representation,
                                                 path_save_fig=str(output_dir / "three_plot_output_ginette.png"), lim_depth=abs(lim_depth),barre_vertical = [])
    if visualisation_propriete_geophy_2D:
        fct.three_plot_propriete_geophy_dharrma(dossier_actuel,debut_representation,fin_representation,pas_representation, facies,
                                                 barre_vertical=[25.5,35.4],lim_depth=abs(lim_depth),path_save_fig =str(output_dir / "three_plot_propriete_geophy.png"))
    if visualisation_observable_geophy_2D:
        fct.three_plot_observable_geophy_dharrma(dossier_actuel,debut_representation,fin_representation,pas_representation, facies, 
                                                 barre_vertical=[25.5,35.4],lim_depth=abs(lim_depth),path_save_fig =str(output_dir / "three_plot_observable_geophy.png"))
    if visualisation_observable_geophy_profil:
        if not required_vp_vs_file.exists():
            print(
                f"Fichier manquant: {required_vp_vs_file} -> relancer avec sismic=True "
                "ou desactiver visualisation_observable_geophy_profil"
            )
        elif homogeneite:
            fct.plot_profil_observable_dharrma(dossier_actuel,jour_profil,facies,representation = representation ,path_fig = str(output_dir / "plot_profil_observable.png"))
        else:
            fct.plot_profil_observable_dharrma(dossier_actuel,jour_profil,facies1,representation = representation ,path_fig = str(output_dir / "plot_profil_observable.png"),facies2 = facies2)
    if visualisation_propriete_geophy_profil:
        if not required_vp_vs_file.exists():
            print(
                f"Fichier manquant: {required_vp_vs_file} -> relancer avec sismic=True "
                "ou desactiver visualisation_propriete_geophy_profil"
            )
        elif homogeneite:
            fct.plot_profil_propriete_dharrma(dossier_actuel,jour_profil,facies,representation = representation,path_fig = str(output_dir / "plot_profil_propriete.png"),facies2 = None)
        else:
            fct.plot_profil_propriete_dharrma(dossier_actuel,jour_profil,facies1,representation = representation,path_fig = str(output_dir / "plot_profil_propriete.png"),facies2 = facies2)
    if visualisation_propriete_hydro_profil:
        fct.plot_sat_jours(dossier_actuel,jour_profil,lim_depth = -2.1)

    if visualisation_wt:
        fct.plot_wt_time(dossier_actuel,pas_jour_x=10)

    if comparaison_wt_piezo:
        fct.plot_comparaison_wt_piezo(
            dossier_actuel,
            path_piezo,
            cote_ngf_piezo,
            date_debut_simulation,
            nbr_jour,
            pas_jour_x=10,
        )

    # Agg backend: figures are generated and saved to disk by plotting functions
    # They are not displayed interactively but are available for review