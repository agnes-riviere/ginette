#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script de préparation des conditions initiales et aux limites pour une simulation Ginette (cas réel)
===============================================================================================
22/07/2026
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
- Toute modification des données d'observation ou des paramètres (pas de temps, durée, etc.) doit être faite dans le bloc PARAMÈTRES UTILISATEUR tout en haut de ce fichier, avant de lancer la simulation.
- Ce script garantit la reproductibilité et la traçabilité de la préparation des cas réels.
- Compatible Windows/Mac/Linux : tous les chemins sont construits avec pathlib, aucun chemin n'est codé en dur.

@author: Agnès Rivière, Samuel Larance, Alexandrine Gesret
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

# Bascule automatique vers un backend non interactif (Agg) si aucun affichage
# n'est disponible (job cluster/SLURM lancé en batch, sans X11) - sinon
# plt.show() plante avec "no display name and no $DISPLAY environment
# variable" au lieu de simplement ne rien afficher. Restreint à Linux : sur
# Mac/Windows l'absence de DISPLAY ne veut pas dire l'absence d'écran (backend
# natif différent), donc on ne veut pas désactiver l'affichage local à tort.
import matplotlib
if not os.environ.get("MPLBACKEND") and not os.environ.get("DISPLAY") and sys.platform.startswith("linux"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Import your modules directly from src_python
sys.path.insert(0, str(project_root / "src" / "src_python"))

from src.src_python import Read_obs

# =============================================================================
# PARAMÈTRES UTILISATEUR — à modifier selon le cas d'étude
# =============================================================================
# Point d'observation (doit correspondre à un dossier dans OBS_point/, ex:
# "lomos230", "Point3_540_SOULTZ", "Point1Touques") : partagé via
# config_lomos.py avec 2_run_real_case.py/3_misfit.py/4_plot_results.py, pour
# que results/{POINT_NAME}/ reste cohérent partout et qu'un autre point ne
# vienne pas écraser les résultats de celui-ci.
sys.path.insert(0, str(Path(__file__).resolve().parent))
# Date de début, durée et pas de temps de la simulation vivent dans
# config_lomos.py (partagé avec 2_run_real_case.py) - période de calage
# recommandée par stallman_diffusivity.py (42j, voir mémoire projet).
from config_lomos import (POINT_NAME, DATE_SIMUL_BG as date_simul_bg,
                           NB_DAY as nb_day, DT as dt)

# État de la simulation : 0 = permanent, 1 = transitoire
state = 1

# Coefficients de traitement des données de pression
# deltaP (H_corrige_24h) pour lomos230 est enregistré en CENTIMÈTRES, pas en mètres :
# coef=0.01 convertit vers les mètres attendus par initial_conditions/boundary_conditions.
coef = 0.01  # coefficient d'échelle pour la pression (cm -> m pour lomos230)
offset = 0  # décalage (correction de référence)

# TempMolo est en eau libre 1-2cm au-dessus du lit, pas dans le sédiment -
# ce n'est pas la même grandeur physique que la T° du milieu poreux à z=0.
# CL haute/basse = TOP_SENSOR/BOTTOM_SENSOR (config_lomos.py), à choisir
# librement parmi SENSOR_DEPTHS. Deux presets déjà étudiés sur lomos230 :
# - Config A (TOP_SENSOR='TempMolo', BOTTOM_SENSOR='Temp4') : TempMolo BRUT en
#   CL haute, domaine naturel 0 à -40cm, calage sur Temp1/Temp2/Temp3.
# - Config C (TOP_SENSOR='Temp1', BOTTOM_SENSOR='Temp4') : Temp1 ET Temp4 tels
#   quels en CL (vraies mesures, aucune correction), calage sur Temp2/Temp3
#   seulement, deltaP ramené au sous-domaine 10-40cm. Résultat le plus proche
#   de l'analyse Stallman indépendante.
# Tout vit dans config_lomos.py, partagé avec 2_run_real_case.py et
# 3_misfit.py, pour ne plus avoir à le resynchroniser à la main partout.
from config_lomos import (TOP_SENSOR, BOTTOM_SENSOR, COMPARISON_SENSORS, SENSOR_DEPTHS,
                           DELTAP_SENSOR, Z_TOP, Z_BOTTOM, DOMAIN_LENGTH, DZ_OBS)

# Géométrie du domaine 1D vertical (colonne), dérivée de TOP_SENSOR/BOTTOM_SENSOR
z_top, z_bottom = Z_TOP, Z_BOTTOM
dz = 0.01        # taille de maille verticale [m] (résolution fine pour les gradients thermiques)
dz_obs = DZ_OBS  # espacement réel entre capteurs [m] (dérivé de SENSOR_DEPTHS)

az = abs(z_top - z_bottom)  # hauteur totale de la colonne [m] (calculée automatiquement)
# =============================================================================


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
    # (fichiers écrits ligne à ligne dans le même ordre que obs_temp : on réutilise donc
    # son index de dates pour éviter que matplotlib n'affiche un axe temporel démarrant en 1970)
    temp_bc = pd.read_csv(path_temp_t, sep=' ', header=None, names=['T_top', 'T_bottom'])
    charge_bc = pd.read_csv(path_charge_t, sep=' ', header=None, names=['h_top', 'h_bottom'])
    temp_bc.index = obs_temp.index[:len(temp_bc)]
    charge_bc.index = obs_temp.index[:len(charge_bc)]

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

    axs[1].set_xlabel('Date')
    plt.tight_layout()
    plt.show()


# =============================================================================
# CHEMINS ET RÉPERTOIRES (détection automatique, compatible Windows/Mac/Linux)
# =============================================================================
importlib.reload(Read_obs)

# Recherche les répertoires src/ et application/ du projet à partir de ce fichier
# (aucun chemin absolu codé en dur, fonctionne quel que soit l'OS ou l'emplacement du dépôt)
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
            src_py = cand_src
        if cand_src2.exists() and src_py is None:
            src_py = cand_src2
        if cand_app.exists():
            app_dir = cand_app
        if src_py or app_dir:
            repo_root = p
            break
        p = p.parent
    return repo_root, src_py, app_dir

REPO_ROOT, SRC_PY, BASE_APP_DIR = find_project_paths()

# Repli sur des chemins relatifs raisonnables si rien n'a été trouvé
if REPO_ROOT is None:
    REPO_ROOT = Path(__file__).resolve().parents[2]
if SRC_PY is None:
    SRC_PY = REPO_ROOT / "src" / "src_python"
if BASE_APP_DIR is None:
    BASE_APP_DIR = REPO_ROOT / "application" / "1D_Stream_aquifer_GridSearch"

# Ajoute src_python au sys.path si disponible
if SRC_PY.is_dir() and str(SRC_PY) not in sys.path:
    sys.path.insert(0, str(SRC_PY))

GINETTE_SENSI = BASE_APP_DIR / "GINETTE_SENSI"
# Un sous-dossier par point (results/{POINT_NAME}/) : un autre point lancé
# ensuite n'écrase pas les résultats de celui-ci.
RESULTS_DIR = BASE_APP_DIR / "results" / POINT_NAME
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Dossier des données d'observation (température, pression) du point choisi ci-dessus
Obs_data = BASE_APP_DIR / "OBS_point" / POINT_NAME

# debug listing to verify (non-blocking)
if SRC_PY.is_dir():
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
                               write_coordonnee_file,
                               run_direct_model,
                               smooth_square_wave)
    from Read_obs import (process_obs_data, convert_dates,read_csv_with_multiple_separators,read_lomos_csv)
    from Plot import plot_initial_conditions
except Exception:
     # fallback to legacy module name if present
    from Direct_model import (setup_ginette2,
                                       initial_conditions,
                                       boundary_conditions,
                                       write_coordonnee_file,
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
for file in Obs_data.glob('*.csv'):
    print(f"- {file}")


# === TRAITEMENT DES DONNÉES D'OBSERVATION ===
# On lit les données de température et pression à partir des fichiers CSV du dossier d'observations.
# Le fichier utilisé est celui du point choisi via POINT_NAME (bloc PARAMÈTRES UTILISATEUR en haut du fichier).
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

# === CL HAUTE/BASSE = TOP_SENSOR/BOTTOM_SENSOR (config_lomos.py) ===
# Toujours des vraies mesures, aucune correction/extrapolation. deltaP est
# mesuré entre la rivière et DELTAP_SENSOR (profondeur fixe de l'installation) -
# si le domaine simulé est plus court, on ne garde que la fraction du gradient
# qui tombe dedans, gradient hydraulique supposé homogène sur toute la colonne
# (pas d'autre mesure dispo).
_needed_cols = set(COMPARISON_SENSORS) | {TOP_SENSOR, BOTTOM_SENSOR, 'deltaP'}
if not _needed_cols.issubset(obs_temp.columns):
    raise ValueError(f"Colonnes manquantes dans les données de {POINT_NAME} pour "
                      f"TOP_SENSOR={TOP_SENSOR!r}/BOTTOM_SENSOR={BOTTOM_SENSOR!r} : "
                      f"{_needed_cols - set(obs_temp.columns)}")

obs_temp['T_top'] = obs_temp[TOP_SENSOR]
obs_temp['T_bottom'] = obs_temp[BOTTOM_SENSOR]
# boundary_conditions() (Direct_model.py) teste 'deltaP' EN PREMIER, avant
# 'h_top'/'h_bottom' - écrire h_top/h_bottom seuls ne suffit pas, ça reste
# ignoré tant que la colonne 'deltaP' existe. On écrase donc deltaP
# directement (bug trouvé le 2026-07-24 : la config C a tourné toute une
# grille avec le deltaP brut, pas le x0.75, avant ce fix).
_deltap_scale = DOMAIN_LENGTH / SENSOR_DEPTHS[DELTAP_SENSOR]
obs_temp['deltaP'] = obs_temp['deltaP'] * _deltap_scale
obs_temp['h_top'] = obs_temp['deltaP']
obs_temp['h_bottom'] = 0.0
# maille1/2/3 de Ginette tombent sur les vraies profondeurs de COMPARISON_SENSORS
# (TOP_SENSOR est consommé comme CL) - il faut renommer AVANT d'écrire
# observed_data.txt sinon 3_misfit.py compare des profondeurs différentes
# entre simulé et observé. Identité si TOP_SENSOR='TempMolo' (Config A).
_orig = obs_temp[COMPARISON_SENSORS].copy()
for _i, _sensor in enumerate(COMPARISON_SENSORS, start=1):
    obs_temp[f'Temp{_i}'] = _orig[_sensor]

print(f"\nCL haute = {TOP_SENSOR}, CL basse = {BOTTOM_SENSOR} (mesures brutes).")
print(f"Calage sur {COMPARISON_SENSORS} (renommés Temp1/Temp2/Temp3 dans observed_data.txt).")
if _deltap_scale != 1.0:
    print(f"deltaP ramené au sous-domaine ({DOMAIN_LENGTH*100:.0f}cm, x{_deltap_scale:.2f}, "
          "gradient homogène supposé).")

# On prépare un DataFrame simplifié pour sauvegarder les données traitées (utile pour vérification ou post-traitement)
data = pd.DataFrame()
data['Time'] = time_vector
data['Temp1'] = obs_temp['Temp1'].to_numpy()
data['Temp2'] = obs_temp['Temp2'].to_numpy()
data['Temp3'] = obs_temp['Temp3'].to_numpy()
data['dates'] = obs_temp['dates'].to_numpy()
# Sauvegarde dans le dossier results
data.to_csv(RESULTS_DIR / "observed_data.txt", sep=" ")

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
#
# NB : initial_conditions() ne regarde qu'un seul instant (all_data.iloc[0]), donc le
# profil peut être localement non monotone avec la profondeur (les capteurs profonds
# sont déphasés par rapport à la surface). On NE lisse PAS cette valeur : l'objectif est
# de comparer le modèle aux observations réelles, pas à une version filtrée. Le modèle
# transitoire (conduction + advection) efface cet écart en environ un temps de diffusion
# sur la colonne (~1 semaine ici) : exclure ce spin-up initial de la comparaison
# modèle/observations plutôt que de modifier l'état initial lui-même.
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




# E_coordonnee.dat n'est normalement écrit que par generate_zone_parameters()
# (2_run_real_case.py) : régénéré ici à la taille du domaine courant
# (z_bottom, dz) avant le plot, sinon plot_initial_conditions() peut lire un
# fichier resté à la taille d'un run précédent (ex: Config A 40 mailles vs
# Config C 30) et planter sur un mismatch de forme entre temp_init et coord.
write_coordonnee_file(z_bottom, dz)

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