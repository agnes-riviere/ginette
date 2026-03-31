#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script de préparation des conditions initiales et aux limites pour une simulation Ginette (cas réel)
===============================================================================================

Ce script automatise la création de tous les fichiers d'entrée nécessaires à une simulation hydrogéologique réaliste avec le code Ginette.
Il s'appuie sur des données de terrain (température, pression) enregistrées dans un fichier CSV d'observations.

**Résumé de la démarche :**
1. **Lecture des données d'observation**
    - Le script lit les mesures de température et pression dans le dossier OBS_point (ex : Point3_540_SOULTZ.csv).
    - Les données sont interpolées sur le pas de temps du modèle (ex : 15 min).

2. **Traitement et sauvegarde des données**
    - Les données traitées sont sauvegardées pour vérification (results/observed_data.txt).

3. **Création des fichiers d'entrée pour Ginette**
    - Le script se place dans le dossier GINETTE_SENSI.
    - Il crée tous les fichiers nécessaires à la simulation :
      - Maillage, paramètres physiques, etc. (via setup_ginette)
      - Conditions initiales (E_temperature_initiale.dat, E_charge_initiale.dat)
      - Conditions aux limites temporelles (E_temp_t.dat, E_charge_t.dat)
      - Fichier principal des conditions aux limites (E_cdt_aux_limites.dat)

4. **Principe de fonctionnement**
    - Les fonctions initial_conditions et boundary_conditions utilisent directement les données du CSV d'observation pour générer les fichiers d'entrée.
    - Ces fichiers sont ensuite utilisés par le script 2_run_real_case.py pour lancer la simulation Ginette.

**À retenir :**
- Toute modification des données d'observation ou des paramètres (pas de temps, durée, etc.) doit être faite ici avant de lancer la simulation.
- Ce script garantit la reproductibilité et la traçabilité de la préparation des cas réels.

@author: Agnès Rivière, Samuel Larance
"""


# IMPORT:

import sys
from pathlib import Path
# Add project root to path
project_root = Path(__file__).resolve().parents[2]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
import importlib
import os
import numpy as np
import pandas as pd
from time import time
import shutil
import multiprocessing as mp
import matplotlib.pyplot as plt

import glob
# Import your modules directly from src_python
sys.path.insert(0, str(project_root / "src" / "src_python"))

from src.src_python import Read_obs

# =============================================================================
# FONCTION DE VISUALISATION DES CONDITIONS AUX LIMITES ET DES DONNÉES BRUTES
# =============================================================================
def plot_boundary_conditions_and_raw_data(obs_temp, path_temp_t, path_charge_t):
    """
    Affiche les séries temporelles des conditions aux limites (température et charge)
    utilisées dans Ginette, ainsi que les données brutes d'observation d'où elles sont issues.
    - obs_temp : DataFrame des observations (issu du CSV)
    - path_temp_t : chemin vers E_temp_t.dat
    - path_charge_t : chemin vers E_charge_t.dat
    """
    import matplotlib.pyplot as plt
    import pandas as pd

    # Lecture des conditions aux limites générées
    temp_bc = pd.read_csv(path_temp_t, sep=' ', header=None, names=['T_top', 'T_bottom'])
    charge_bc = pd.read_csv(path_charge_t, sep=' ', header=None, names=['h_top', 'h_bottom'])

    # On suppose que obs_temp contient les colonnes 'T_top', 'T_bottom', 'deltaP' ou 'h_top', 'h_bottom'
    # Si ce n'est pas le cas, on affiche toutes les colonnes disponibles

    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # --- Température ---
    if 'T_top' in obs_temp.columns and 'T_bottom' in obs_temp.columns:
        axs[0].plot(obs_temp.index, obs_temp['T_top'], label='T_top (brut)', color='tab:blue', alpha=0.5)
        axs[0].plot(obs_temp.index, obs_temp['T_bottom'], label='T_bottom (brut)', color='tab:orange', alpha=0.5)
    elif 'TempMolo' in obs_temp.columns and 'Temp4' in obs_temp.columns:
        axs[0].plot(obs_temp.index, obs_temp['TempMolo'], label='TempMolo (brut)', color='tab:blue', alpha=0.5)
        axs[0].plot(obs_temp.index, obs_temp['Temp4'], label='Temp4 (brut)', color='tab:orange', alpha=0.5)
    else:
        for col in obs_temp.columns:
            if 'temp' in col.lower():
                axs[0].plot(obs_temp.index, obs_temp[col], label=f'{col} (brut)', alpha=0.5)

    axs[0].plot(temp_bc.index, temp_bc['T_top'], label='T_top (BC)', color='tab:blue', linestyle='--')
    axs[0].plot(temp_bc.index, temp_bc['T_bottom'], label='T_bottom (BC)', color='tab:orange', linestyle='--')
    axs[0].set_ylabel('Température (°C)')
    axs[0].set_title('Conditions aux limites de température et données brutes')
    axs[0].legend()

    # --- Charge/pression ---
    if 'deltaP' in obs_temp.columns:
        axs[1].plot(obs_temp.index, obs_temp['deltaP'], label='deltaP (brut)', color='tab:green', alpha=0.5)
    elif 'h_top' in obs_temp.columns and 'h_bottom' in obs_temp.columns:
        axs[1].plot(obs_temp.index, obs_temp['h_top'], label='h_top (brut)', color='tab:green', alpha=0.5)
        axs[1].plot(obs_temp.index, obs_temp['h_bottom'], label='h_bottom (brut)', color='tab:red', alpha=0.5)
    else:
        for col in obs_temp.columns:
            if 'h' in col.lower() or 'press' in col.lower():
                axs[1].plot(obs_temp.index, obs_temp[col], label=f'{col} (brut)', alpha=0.5)

    axs[1].plot(charge_bc.index, charge_bc['h_top'], label='h_top (BC)', color='tab:green', linestyle='--')
    axs[1].plot(charge_bc.index, charge_bc['h_bottom'], label='h_bottom (BC)', color='tab:red', linestyle='--')
    axs[1].set_ylabel('Charge/Pression')
    axs[1].set_title('Conditions aux limites de charge et données brutes')
    axs[1].legend()

    axs[1].set_xlabel('Pas de temps (index)')
    plt.tight_layout()
    plt.show()


# Set up directories and data paths

# Get current script directory
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = str(SCRIPT_DIR)

# Navigate to sibling directories
GINETTE_SENSI = os.path.join(BASE_DIR, "GINETTE_SENSI")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
# This path contains the observational data (temperature sensors, pressure measurements)
Obs_data = os.path.join(BASE_DIR,'OBS_point/Point3_540_SOULTZ/')
# =============================================================================
# SIMULATION SETUP: Define temporal parameters and read observational data
# =============================================================================
importlib.reload(Read_obs)
# Start date and time for the simulation
# This corresponds to when field measurements began
date_simul_bg = pd.to_datetime("2025/09/26 12:00:00")

# Simulation state configuration:
# 0 = steady state (time-independent, equilibrium conditions)
# 1 = transient state (time-dependent, dynamic evolution)
# We use transient state to capture temporal variations in temperature and flow
state = 1

# Data processing coefficients
coef = 1    # Scaling coefficient for pressure measurements
offset = 0  # Offset for pressure measurements (baseline correction)
# Simulation duration in days
# 30 days provides sufficient time to observe seasonal patterns and model spin-up
nb_day = 30
# =============================================================================
# MODEL GEOMETRY AND DISCRETIZATION SETUP
# =============================================================================

# TEMPORAL DISCRETIZATION
# Time step in seconds (900s = 15 minutes)
# This matches the measurement frequency and ensures numerical stability
dt = 900

# SPATIAL DOMAIN DEFINITION (1D vertical column)
z_top = 0.0      # Surface elevation (stream bed) [m]
z_bottom = -0.4  # Bottom of model domain [m] (40 cm below streambed)
az = abs(z_top - z_bottom)  # Total column height [m]

# GRID DISCRETIZATION
dz = 0.01        # Vertical cell size [m] (1 cm resolution)
                 # Fine discretization needed to capture thermal gradients

# OBSERVATION DEPTHS
dz_obs = 0.1     # Spacing between temperature sensors [m] (10 cm)
                 # Sensors at -10, -20, -30, -40 cm depths


# =============================================================================


# Find project-relative src and application directories (no absolute paths)
def find_project_paths(start_file=__file__):
    p = Path(start_file).resolve().parent
    repo_root = None
    src_py = None
    app_dir = None
    for _ in range(8):
        # candidate locations relative to current parent
        cand_src = p / "src" / "src_python"
        cand_src2 = p / "src_python"
        cand_app = p / "application" / "1D_Stream_aquifer_GridSearch"
        if cand_src.exists():
            src_py = str(cand_src)
        if cand_src2.exists() and src_py is None:
            src_py = str(cand_src2)
        if cand_app.exists():
            app_dir = str(cand_app)
        if src_py or app_dir:
            repo_root = str(p)
            break
        p = p.parent
    return repo_root, src_py, app_dir

REPO_ROOT, SRC_PY, BASE_APP_DIR = find_project_paths()

# Fallback to reasonable relative defaults if not found
if REPO_ROOT is None:
    REPO_ROOT = Path(__file__).resolve().parents[2].as_posix()
if SRC_PY is None:
    SRC_PY = os.path.join(REPO_ROOT, "src", "src_python")
if BASE_APP_DIR is None:
    BASE_APP_DIR = os.path.join(REPO_ROOT, "application", "1D_Stream_aquifer_GridSearch")

# add src python to sys.path if available
if os.path.isdir(SRC_PY) and SRC_PY not in sys.path:
    sys.path.insert(0, SRC_PY)
if BASE_APP_DIR is None:
    BASE_APP_DIR = os.path.join(REPO_ROOT, "application", "1D_Stream_aquifer_GridSearch")

# --- Ajout: dossier de sortie local à l'application (modifiable) ---
RESULTS_DIR = os.path.join(BASE_APP_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# debug listing to verify (non-blocking)
if os.path.isdir(SRC_PY):
    try:
        print(f"Using src_python: {SRC_PY}")
        print("Files:", sorted(os.listdir(SRC_PY)))
    except Exception:
        pass

# Import project modules from src_python (robust to module name/case)
try:
    # preferred: modules as they appear in src/src_python
    from Direct_model import (setup_ginette2,setup_ginette,
                               initial_conditions,
                               boundary_conditions,
                               run_direct_model,
                               smooth_square_wave)
    from Read_obs import (process_obs_data, convert_dates,read_csv_with_multiple_separators)
    from Plot import plot_initial_conditions
except Exception:
     # fallback to legacy module name if present
    from Direct_model import (setup_ginette2,
                                       initial_conditions,
                                       boundary_conditions,
                                       run_direct_model,
                                       smooth_square_wave)


try:
    from Init_folders import prepare_ginette_directories
except Exception:
    # fallback if an alternative utils module exists
    try:
        from utils_ginette import prepare_ginette_directories
    except Exception:
        # minimal fallback: create directory helper
        def prepare_ginette_directories(path):
            os.makedirs(path, exist_ok=True)

# provide a portable copy_file helper if the project does not expose one
try:
    from Init_folders import copy_file  # some versions may provide it
except Exception:
    def copy_file(src, dst_dir):
        """Copy src into dst_dir (create dst_dir if needed)."""
        os.makedirs(dst_dir, exist_ok=True)
        if not os.path.exists(src):
            raise FileNotFoundError(f"Source not found: {src}")
        shutil.copy(src, os.path.join(dst_dir, os.path.basename(src)))







# List available CSV files in the observational data directory
print("CSV files in observational data directory:")
for file in glob.glob(os.path.join(Obs_data, '*.csv')):
    print(f"- {file}")


# === TRAITEMENT DES DONNÉES D'OBSERVATION ===
# On lit les données de température et pression à partir des fichiers CSV du dossier d'observations.
# Le fichier principal utilisé est : OBS_point/Point3_540_SOULTZ/Point3_540_SOULTZ.csv
# La fonction process_obs_data :
#   - lit les séries temporelles de température et pression,
#   - interpole les données sur le pas de temps du modèle (ici 15 min),
#   - applique un contrôle qualité,
#   - retourne un DataFrame synchronisé prêt à être utilisé pour les conditions initiales et aux limites.
obs_temp = process_obs_data(Obs_data, date_simul_bg, coef, offset, nb_day)
# (Remarque : pour changer le pas de temps, modifier la variable dt plus haut)
print(f"\nObservational data loaded successfully:")
print(f"- Time period: {obs_temp.index.min()} to {obs_temp.index.max()}")
print(f"- Number of time steps: {len(obs_temp)}")
print(f"- Available measurements: {list(obs_temp.columns)}")
print(f"- Data shape: {obs_temp.shape}")

# On ajoute une colonne 'Time' (secondes depuis le début de la simulation)
time_vector = (obs_temp.index - date_simul_bg).total_seconds().to_numpy()
obs_temp.insert(0, 'Time', time_vector)

# On prépare un DataFrame simplifié pour sauvegarder les données traitées (utile pour vérification ou post-traitement)
data = pd.DataFrame()
data['Time'] = time_vector
data['Temp1'] = obs_temp['Temp1'].to_numpy()
data['Temp2'] = obs_temp['Temp2'].to_numpy()
data['Temp3'] = obs_temp['Temp3'].to_numpy()
data['dates'] = obs_temp['dates'].to_numpy()
# Sauvegarde dans le dossier results
data.to_csv(os.path.join(RESULTS_DIR, "observed_data.txt"), sep=" ")

# === CRÉATION DES FICHIERS D'ENTRÉE POUR GINETTE ===
# On se place dans le dossier GINETTE_SENSI où seront créés les fichiers d'entrée du modèle Ginette.
os.chdir(GINETTE_SENSI)

# La fonction setup_ginette crée les fichiers de paramètres principaux du modèle (maillage, paramètres physiques, etc.)
# Elle retourne aussi la liste des profondeurs d'observation utilisées.
z_obs = setup_ginette(dt, state, nb_day, z_top, z_bottom, az, dz, date_simul_bg, dz_obs)

print(f"Model domain setup:")
print(f"- Vertical extent: {z_top} to {z_bottom} m")
print(f"- Cell size: {dz} m")
print(f"- Number of cells: {int(az/dz)}")
print(f"- Time step: {dt} s ({dt/60} minutes)")
print(f"- Observation depths: {z_obs} m")



# =============================================================================
# INITIAL AND BOUNDARY CONDITIONS SETUP
# =============================================================================


# === CONDITIONS INITIALES ===
# Cette étape crée les fichiers d'initialisation du modèle Ginette à t=0, à partir des observations :
#   - E_temperature_initiale.dat : profil de température initial interpolé sur la colonne
#   - E_charge_initiale.dat      : profil de charge initiale (pression) interpolé
# Les données proviennent du DataFrame obs_temp, issu du CSV d'observations.
initial_conditions(obs_temp, z_top, z_bottom, dz, z_obs)
# (Voir src/src_python/Direct_model.py pour le détail de la fonction)

# === CONDITIONS AUX LIMITES ===
# Cette étape crée les fichiers de conditions aux limites (hydrauliques et thermiques) pour la simulation transitoire :
#   - E_cdt_aux_limites.dat : fichier principal des conditions aux limites (format Ginette)
#   - E_charge_t.dat        : conditions aux limites de charge (pression) en fonction du temps
#   - E_temp_t.dat          : conditions aux limites de température en fonction du temps
# Ces fichiers sont générés à partir des observations (CSV) et du pas de temps choisi.
boundary_conditions(obs_temp, dt)
# (Voir src/src_python/Direct_model.py pour le détail de la fonction)

print("Initial and boundary conditions applied:")
print("- Initial temperature profile interpolated from observations")
print("- Initial pressure profile derived from surface measurements")
print("- Time-varying surface conditions from observational data")


# =============================================================================
# PLOT INITIAL CONDITIONS
# =============================================================================




# Plot the initial temperature and pressure profiles
# This visualization helps verify that the initial conditions are realistic
# and properly interpolated from the observational data
plot_initial_conditions()
# Affiche la figure à l'écran (important si le script est lancé en mode batch ou IDE)

# The plots show:
# 1. Initial temperature profile: Interpolated from sensor measurements
# 2. Initial pressure profile: Derived from surface pressure measurements
# 
# These initial conditions serve as the starting point (t=0) for the simulation
# and should represent a reasonable approximation of the actual field conditions



plot_boundary_conditions_and_raw_data(
        obs_temp,
        "E_temp_t.dat",
        "E_charge_t.dat")
plt.show()