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
- Toute modification des données d'observation ou des paramètres (pas de temps, durée, etc.) doit être faite dans le bloc PARAMÈTRES UTILISATEUR tout en haut de ce fichier, avant de lancer la simulation.
- Ce script garantit la reproductibilité et la traçabilité de la préparation des cas réels.
- Compatible Windows/Mac/Linux : tous les chemins sont construits avec pathlib, aucun chemin n'est codé en dur.

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
# Nom du point d'observation : doit correspondre à un dossier dans OBS_point/
# (ex: "lomos230", "Point3_540_SOULTZ", "Point1Touques")
POINT_NAME = "lomos230"

# Date de début de la simulation (doit être dans la plage de données du point choisi)
# Période de calage recommandée par l'analyse Stallman/Cheviron (voir
# 0_stallman_diffusivity_lomos230.py) : la plus longue fenêtre continue où la
# méthode thermique et le gradient de charge s'accordent sur le sens d'écoulement.
date_simul_bg = pd.to_datetime("2026/05/15 04:00:00")

# État de la simulation : 0 = permanent, 1 = transitoire
state = 1

# Coefficients de traitement des données de pression
# deltaP (H_corrige_24h) pour lomos230 est enregistré en CENTIMÈTRES, pas en mètres :
# coef=0.01 convertit vers les mètres attendus par initial_conditions/boundary_conditions.
coef = 0.01  # coefficient d'échelle pour la pression (cm -> m pour lomos230)
offset = 0  # décalage (correction de référence)

# Durée de simulation en jours
nb_day = 42

# Pas de temps en secondes (900 s = 15 min, doit correspondre au pas des observations)
dt = 900

# TempMolo (capteur de "surface") est mesuré 1-2cm au-dessus du lit, pas au
# contact du sédiment, ce qui sur-amplifie le cycle diurne imposé comme
# condition de Dirichlet à z=0 pour lomos230/231. USE_TEMP1_AS_TOP_BC corrige
# ça en décalant tout le schéma d'un cran de profondeur (Temp1 mesuré devient
# la condition limite haute). Vérifier avec l'analyse Stallman avant de
# l'activer pour un nouveau point : voir OBS_point/notes_points_obs.txt pour
# le diagnostic complet, les alternatives écartées et les résultats de misfit
# par point.
USE_TEMP1_AS_TOP_BC = POINT_NAME in ("lomos230", "lomos231")
CORRECT_TEMPMOLO_SURFACE = False
DZ_ANCHOR_TO_SURFACE = 0.10  # distance (m) Temp1 -> surface du sédiment (z=0)

# Géométrie du domaine 1D vertical (colonne)
# Avec USE_TEMP1_AS_TOP_BC, le domaine est décalé de -10cm : le "haut" du
# modèle correspond à la profondeur réelle de Temp1, pas à la surface du lit.
_DEPTH_SHIFT = 0.10 if USE_TEMP1_AS_TOP_BC else 0.0
z_top = 0.0 - _DEPTH_SHIFT       # altitude de la surface (lit de rivière) [m]
z_bottom = -0.4 - _DEPTH_SHIFT   # bas du domaine [m]
dz = 0.01        # taille de maille verticale [m] (résolution fine pour les gradients thermiques)
dz_obs = 0.1     # espacement entre les capteurs de température [m]

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
RESULTS_DIR = BASE_APP_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

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
                               run_direct_model,
                               smooth_square_wave)
    from Read_obs import (process_obs_data, convert_dates,read_csv_with_multiple_separators,read_lomos_csv)
    from Plot import plot_initial_conditions
    from Stallman_analysis import (fit_sinusoid, fit_local_attenuation_phase,
                                    extrapolate_attenuation_phase, predict_true_surface_signal,
                                    predict_deeper_signal)
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

# === DÉCALAGE : TEMP1 COMME CONDITION LIMITE HAUTE (à la place de TempMolo) ===
# Voir USE_TEMP1_AS_TOP_BC dans PARAMÈTRES UTILISATEUR et
# OBS_point/notes_points_obs.txt pour le diagnostic complet. On décale toutes
# les colonnes de température d'un cran de profondeur : Temp1 (mesure réelle)
# devient la nouvelle "TempMolo" (condition de Dirichlet du haut), Temp2->Temp1,
# Temp3->Temp2, Temp4->Temp3 (nouveaux points de comparaison, tous mesurés).
# La nouvelle condition limite du bas (-50cm, 10cm sous Temp4) est prédite par
# extrapolation physique (predict_deeper_signal) plutôt que réutilisée telle
# quelle depuis Temp4 - la réutilisation brute rendait Temp3 et Temp4
# rigoureusement identiques, un biais structurel (détail dans le fichier de
# notes ci-dessus).
if USE_TEMP1_AS_TOP_BC and all(c in obs_temp.columns for c in ['TempMolo', 'Temp1', 'Temp2', 'Temp3', 'Temp4']):
    _orig = obs_temp[['Temp1', 'Temp2', 'Temp3', 'Temp4']].copy()
    _depths_orig = {'Temp1': 0.10, 'Temp2': 0.20, 'Temp3': 0.30, 'Temp4': 0.40}
    _fits_orig = {name: fit_sinusoid(time_vector, _orig[name].to_numpy(), 86400.0) for name in _depths_orig}
    _local_deep = fit_local_attenuation_phase(_fits_orig, _depths_orig, 86400.0)
    _z_target_deep = _depths_orig['Temp4'] + DZ_ANCHOR_TO_SURFACE / 2  # milieu Temp4 <-> nouvelle frontière
    _a_deep, _b_deep = extrapolate_attenuation_phase(_local_deep, _z_target_deep)
    _deep_pred = predict_deeper_signal(time_vector, _orig['Temp4'].to_numpy(), _fits_orig['Temp4'], 86400.0,
                                        dz_deeper=DZ_ANCHOR_TO_SURFACE, a_pred=_a_deep, b_pred=_b_deep)

    obs_temp['TempMolo_raw'] = obs_temp['TempMolo']  # conservé pour référence, non utilisé
    obs_temp['TempMolo'] = _orig['Temp1']
    obs_temp['Temp1'] = _orig['Temp2']
    obs_temp['Temp2'] = _orig['Temp3']
    obs_temp['Temp3'] = _orig['Temp4']
    obs_temp['Temp4_raw'] = _orig['Temp4']  # conservé pour référence/diagnostic
    obs_temp['Temp4'] = _deep_pred['predicted']
    print("\n=== USE_TEMP1_AS_TOP_BC actif ===")
    print("Condition limite haute = Temp1 mesuré (plus d'extrapolation/correction).")
    print("Points de comparaison (nouveaux Temp1/2/3) = anciens Temp2/Temp3/Temp4, tous mesurés.")
    print(f"Condition limite basse (-50cm) prédite par extrapolation depuis Temp4 : "
          f"amplitude diurne {_deep_pred['amplitude_raw']:.4f} -> {_deep_pred['amplitude_pred']:.4f} degC")

# === CORRECTION DU SIGNAL DE SURFACE (TempMolo) ===
# Voir le commentaire de CORRECT_TEMPMOLO_SURFACE dans PARAMÈTRES UTILISATEUR.
if CORRECT_TEMPMOLO_SURFACE and 'TempMolo' in obs_temp.columns and 'Temp1' in obs_temp.columns:
    depths_interior = {'Temp1': 0.10, 'Temp2': 0.20, 'Temp3': 0.30, 'Temp4': 0.40}
    fits_interior = {name: fit_sinusoid(time_vector, obs_temp[name].to_numpy(), 86400.0)
                     for name in depths_interior if name in obs_temp.columns}

    local_att_phase = fit_local_attenuation_phase(fits_interior, depths_interior, 86400.0)
    z_target = depths_interior['Temp1'] - DZ_ANCHOR_TO_SURFACE / 2  # milieu de la paire surface-Temp1
    a_pred, b_pred = extrapolate_attenuation_phase(local_att_phase, z_target)

    correction = predict_true_surface_signal(
        time_vector, obs_temp['TempMolo'].to_numpy(), fits_interior['Temp1'], 86400.0,
        dz_to_surface=DZ_ANCHOR_TO_SURFACE, a_pred=a_pred, b_pred=b_pred)

    print("\n=== Correction du signal de surface (TempMolo) ===")
    print(local_att_phase.to_string(index=False))
    print(f"a_pred (extrapolé à z={z_target:.2f}m): {a_pred:.2f} /m")
    print(f"b_pred (extrapolé à z={z_target:.2f}m): {b_pred:.2f} rad/m")
    print(f"Amplitude diurne TempMolo brute : {correction['amplitude_raw']:.4f} degC")
    print(f"Amplitude diurne surface prédite: {correction['amplitude_pred']:.4f} degC "
          f"(ratio {correction['amplitude_pred']/correction['amplitude_raw']:.3f})")

    obs_temp['TempMolo_raw'] = obs_temp['TempMolo'].copy()
    obs_temp['TempMolo'] = correction['corrected']

    fig_corr, ax_corr = plt.subplots(figsize=(10, 4))
    ax_corr.plot(obs_temp.index, obs_temp['TempMolo_raw'], label='TempMolo brut (1-2cm au-dessus du lit)',
                 color='tab:blue', alpha=0.6)
    ax_corr.plot(obs_temp.index, obs_temp['TempMolo'], label='TempMolo corrigé (surface du sédiment, z=0)',
                 color='tab:red', linestyle='--')
    ax_corr.set_ylabel('Température (°C)')
    ax_corr.set_title('Correction du signal de surface utilisé comme condition limite')
    ax_corr.legend()
    fig_corr.tight_layout()
    fig_corr.savefig(RESULTS_DIR / "tempmolo_correction.png")

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