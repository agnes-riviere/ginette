'''
Modèle DHARRMA (Direct HydrogeophysicAl Resistivity and Refraction Modeling Application)


Code permettant de lancer le modèle direct transitoire hydrogéophysique développé par N. RADIC, A. RIVIERE, L. BODET, M. GAUTIER, A. GESRET, R. MARTIN en 2025
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


##################### PARTIE 0 : PARAMETRE DU MODÈLE #########################################

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


# Choix des parties à faire tourner dans le code
thermique = True
sismic = False
electrique = False
visualisation = True



######################## PARTIE I : INPUT DU MODÈLE ##########################################

# Paramètres Généraux

nbr_jour = 120 # nombre de jours de simulation au total
facies = 'silt' # Faciès utilisé dans la simulation (cf Carsel and Parish (1986)) 


###### Paramètre de simulation Hydro GINETTE-----------------------------------------------------------------------------------------------------------------------
lancer_ginette = True
# Maillage de GINETTE
depth_top = 0 # (en m) haut du modèle
depth_bottom = -4 # (en m) bas du modèle
hauteur_WT_initial = -2 # (en m) Hauteur de la nappe à l'état initiale
dz = 0.01 # (en m) largeur des mailles du modèle hydro

# Paramètre hydro
pas_hydro = 900 # En seconde, Pas de temps à utiliser pour ma simulation hydro
k1 = 6.94E-14 #1E-11# Perméabilité
soil = selectSoilType(facies) # Paramètre selon le facies et la distribution de Carsel and Parrish
phi_soil = soil[3] # Porosité
Swr_soil = soil[7] # Saturation résiduelle

#Van genuchten parameter
alpha_soil = soil[4]
nvg_soil = soil[5]



if thermique:
    ith=1
else:
    ith=0
Dm.setup_ginette_DHARRMA(pas_hydro, ith, nbr_jour, depth_top, depth_bottom, abs(depth_top-depth_bottom), dz, pas_hydro,hauteur_WT_initial)

# Création de fichier d'entrée pour GINETTE--------------------------------------------------------------------------------------------------------------------
Creation_temp = True
Creation_infiltation = True
Evapo = False

# Création du fichier de température
# Si un fichier de température existe déjà il est possible de le mettre dans le fichier E_temp_t.dat
if Creation_temp :

    bottom_temp = 10 #(°C) Temp stable en bas du modèle
    deg_per_day = 0.2 #0#(°C) Linear temperature increase over time
    daily_fluctuation = 5
    weekly_fluctuation = 3
    temp_offset = 0

    fct.creation_temp(nbr_jour,pas_hydro,bottom_temp,deg_per_day,daily_fluctuation,weekly_fluctuation,temp_offset)

# Création du fichier d'infiltration
# Si un fichier d'infiltration existe déjà il est possible de le mettre dans le fichier E_debit_haut_t.dat
if Creation_infiltation :
    min_infiltration = 0
    max_infiltration = 4E-08#8 # m/s
    duree_transition = 1 #0.125 #  Jours
    duree_max_infiltration= 1 #0.25#100 # Jours
    debut_jour_infiltration = [7,12,17,22,30,87,92,97,102] #[30,102]#[30,87,95,110]#[i for i in range(1,120)]#[5,35,65,95] #[1]#


    if Evapo :
        max_evapo = 7E-8
        duree_transition_evapo = 0.125 # 1 # Jours
        duree_max_evapo = 0.25 # 1 # Jours
        debut_jour_evapo = [i+0.5 for i in range(1,120)]
        fct.creation_infiltration_evapo(nbr_jour,pas_hydro,min_infiltration,max_infiltration,duree_transition,debut_jour_infiltration,duree_max_infiltration,0,max_evapo,duree_transition_evapo,debut_jour_evapo,duree_max_evapo)
    else :
        fct.creation_infiltration(nbr_jour,pas_hydro,min_infiltration,max_infiltration,duree_transition,debut_jour_infiltration,duree_max_infiltration)

###### Paramètre de simulation Thermique GINETTE--------------------------------------------------------------------------------------------------------------------

lambda1 = 2.4
C1 = 897
rho1 = 2400

Dm.generate_zone_parameters_DHARRMA(thermique,depth_bottom,depth_top, dz, k1, phi_soil, alpha_soil, nvg_soil, Swr_soil, lambda1, C1, rho1)

###### Paramètre de simulation Sismique (avec Geopsy) --------------------------------------------------------------------------------------------------------------
######### Parametre simulation #############
debut_sim_sis = 1 # en jours
pas_sim_sis = 1 # en jours
fin_sim_sis = nbr_jour # (compris) en jours 

dz_sis = 0.01 # Discrétisation verticale (en m)

first_arrival_calcul = True

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
soiltypes = [facies]

# Grains/agregate parameters per layer
Ns = [9] # Coordination Number (number of contact per grain) | default = 8
fracs = [0.3] # Fraction of non-slipping grains (helps making the soil less stiff) | default = 0.3
# Four possible RP models:
# kk = 1 # Constant Pe (see the approach of Zyserman et al., 2017)
# kk = 2 # Pe without suction
# kk = 3 # Pe with suction (cf. Solazzi et al. 2021)
kk = 4 # Données de pression externe

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


###### Paramètre de simulation Electrique (avec PyGImly) -----------------------------------------------------------------------------------------------------------
######### Parametre simulation #############
debut_sim_elec = 1 # en jours
pas_sim_elec = 1 # en jours
fin_sim_elec = nbr_jour # (compris) en jours 

elec_static = True

### Resistivité Vrai : Loi Petrophysique
#Parametre Loi d'Archie
a_archie = 1.196 # Facteur de tortuosité ]0.5;1.5] #limon = 1.196, sable = 1.147
rho_water_25 = 75 # Resistivité elec du fluide ici de l'eau à 25°C (en ohm.m)
a_T = 0.02 # Compensateur de Température Hayashi 2004 et valeur de Matthes 1982
m_archie = 1.929 # Exposant de concentration # limon = 1.929, sable = 2.135
n_archie = 2.338 # Exposant de saturation # limon = 2.338, sable = 0.858
# Parametre a,m et n calculé avec les techniques conventionelles.

beta_s = 5.2E-9

Waxman_smits = True
if Waxman_smits:
    # B_WS = Equation dans waxman smith
    wsand = soil[0]
    wclay = soil[1]
    wsilt = soil[2]
    CEC = wclay*19270 # illite Woodruff and Revil 2011
    Q_v = rho1*((1-phi_soil)/phi_soil)*CEC

    print("Q_v =",Q_v)


### Resistivité mesurée : Problème direct

ab2_scale = 'log' # Recommandée
# ab2_scale = 'linear' 
ab2_min = 1.5 # en mètre
ab2_max = 100 # en mètre
ab2_nbr_pt = 100 # Nombre d'écartement des électrodes AB

mn2 = 1.0 # en mètre, écartement des électrodes MN

###### Paramètre Visualisation -------------------------------------------------------------------------------------------------------------------------------------

visualisation_temp = True
visualisation_pluie = True

#### Visualisation 2D ####
debut_representation = 1 # en jours
fin_representation = nbr_jour # en jours
pas_representation = 1 # en jours

lim_depth = 2.1 # (en mètre) Profondeur jusqu'à laquelle on représente le sol
#Sortie ginette 2D
visualisation_output_ginette = True
#Sortie FWD model
visualisation_propriete_geophy_2D = False
visualisation_observable_geophy_2D = False

#Profil jour
jour_profil = [1,30,40,60] #[i for i in range(1,110,30)]#[10,50,110] #[1,30,60,90,119] # [1,5,10,15,20]#Jour où on plot les profils des différents modèles suivant ce qu'on veut.
representation = 1 # Représentation sur un seul graph
# representation = 2 # Représentation sur plusieurs graph
visualisation_propriete_geophy_profil = False
visualisation_observable_geophy_profil = True

DEBUG = False

#################### PARTIE II : MODÈLE HYDRO/THERMIQUE ######################################

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
    subprocess.call(["./ginette"])
    # Revenez au répertoire précédent
    os.chdir('..')
    print('Run model hydro et thermique Ok')

####################### PARTIE III : MODÈLE SISMIQUE #########################################

if sismic:
    depth = depth_top-depth_bottom
    zs = -np.arange(dz, depth + dz, dz) # Depth positions (negative downward) [m]
    thks = np.diff(np.abs(zs)) # thickness vector [m]

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
        
        for prof in range (len(hs)):
            Swes_val = (saturation_profil[prof] - Swr_soil) / (1 - Swr_soil)
        
            Swes.append(Swes_val)

        # Effective Grain Properties (constant with depth)
        mus, ks, rhos, nus = hillsAverage(mu_clay, mu_silt, mu_sand, rho_clay,
                                            rho_silt, rho_sand, k_clay, k_silt,
                                            k_sand, soiltypes)
        
        thicknesses = [depth]
        
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
            path_vp_vs = f'sismique/{soiltypes[0]}/output_SL_kk{kk}_Vp_Vs.dat'
            path_PS_v = f'sismique/{soiltypes[0]}/output_SL_kk{kk}_PS_v_Phase.dat'                

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

            kappa_s = beta_s*Q_v/sat_profil_array
            kappa = (pow(sat_profil_array,n_archie)/(a_archie*pow(phi_soil,-m_archie)))*((1/rho_water_T)+kappa_s)
            # kappa = (pow(sat_profil_array,n_archie)/(a_archie*pow(phi_soil,-m_archie)))*((1/25)+0.5)
            rho_vrai = 1/kappa
            
            # rho_vrai = a_archie*rho_water_T*pow(phi_soil,-m_archie)*pow(sat_profil_array,-n_archie)*((1+rho_water_T*B_WS*Qv_WS)/(1+(rho_water_T*B_WS*Qv_WS/sat_profil_array)))
            # kappa = 5.2E-9


        else : # Loi d'Archie 
            rho_vrai = a_archie*rho_water_T*pow(phi_soil,-m_archie)*pow(sat_profil_array,-n_archie)

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
            path_rho_vrai = f'elec/{soiltypes[0]}/rho_vrai.dat'
            path_rho_app_AB2 = f'elec/{soiltypes[0]}/rho_app_AB2.dat'
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
    path_saturation = '/home/nradic/Documents/ginette/application/model_dharrma/input_ginette/S_saturation_profil_t.dat'
    # path_temp = '/home/nradic/Documents/ginette/application/model_dharrma/save_data/data_article_1/scenario_5_pluies/input_ginette/S_temperature_t.dat'
    path_temp = '/home/nradic/Documents/ginette/application/model_dharrma/input_ginette/S_temperature_t.dat'
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

        path_sat_static = '/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/saturation_profil.dat'
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


        if i == 0 : # Ecriture rho vrai permanent
            if not os.path.exists('/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_vrai_static_temperature.dat'):
                os.makedirs(os.path.dirname('/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_vrai_static_temperature.dat'), exist_ok=True)
            with open(f'/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_vrai_static_temperature.dat','w') as f:
                for z_value,rho_vrai_value in zip(zs,rho_vrai):
                    f.write(f"{temps_ginette} {z_value} {rho_vrai_value}\n")
        else :
             with open(f'/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_vrai_static_temperature.dat','a') as f:
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

        if i == 0 : # Ecriture rho mesuré permanent
            if not os.path.exists('/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_app_AB2_static_temperature.dat'):
                os.makedirs(os.path.dirname('/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_app_AB2_static_temperature.dat'), exist_ok=True)
            with open(f'/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_app_AB2_static_temperature.dat','w') as f:
                for ab2_value,ra_value in zip(ab2,ra):
                    f.write(f"{temps_ginette} {ab2_value} {ra_value}\n")
        else :
             with open(f'/home/nradic/Documents/ginette/application/model_dharrma/save_data/resultat_static/elec/rho_app_AB2_static_temperature.dat','a') as f:
                for ab2_value,ra_value in zip(ab2,ra):
                    f.write(f"{temps_ginette} {ab2_value} {ra_value}\n")

################### PARTIE V : VISUALISATION DES DONNÉES #####################################

if visualisation :
    if visualisation_temp :
        fct.plot_temp_scenario(f'{dossier_actuel}/input_ginette/E_temp_t.dat')
    if visualisation_pluie:
        fct.plot_graphe_pluie(f'{dossier_actuel}/input_ginette/E_debit_haut_t.dat',nbr_jour, pas_hydro,path_pluie_B = None, 
                              liste_mesure = jour_profil, hauteur_WT_A = None, hauteur_WT_B = None, barre_vertical = [],path_fig =None, max_cumul=1E-3)
    if visualisation_output_ginette:
        fct.three_plot_output_ginette_dharrma(dossier_actuel,debut_representation,fin_representation,pas_representation,
                                                 path_save_fig=None, lim_depth=abs(lim_depth),barre_vertical = [])
    if visualisation_propriete_geophy_2D:
        fct.three_plot_propriete_geophy_dharrma(dossier_actuel,debut_representation,fin_representation,pas_representation, facies,
                                                 barre_vertical=[25.5,35.4],lim_depth=abs(lim_depth),path_save_fig =None)
    if visualisation_observable_geophy_2D:
        fct.three_plot_observable_geophy_dharrma(dossier_actuel,debut_representation,fin_representation,pas_representation, facies, 
                                                 barre_vertical=[25.5,35.4],lim_depth=abs(lim_depth),path_save_fig =None)
    if visualisation_observable_geophy_profil:
        fct.plot_profil_observable_dharrma(dossier_actuel,jour_profil,facies,representation = representation ,path_fig = None)
    if visualisation_propriete_geophy_profil:
        fct.plot_profil_propriete_dharrma(dossier_actuel,jour_profil,facies,representation = representation,path_fig = None)
    
    if True :
        fct.plot_only_sat(dossier_actuel,86400*8,86400*11,900*4*12,lim_depth=-3)

    plt.show()