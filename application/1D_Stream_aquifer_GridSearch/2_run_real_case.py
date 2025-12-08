#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 09:51:06 2025

@author: Maxime Gautier
"""


# IMPORT:
import sys
from pathlib import Path
# Add project root to path
project_root = Path(__file__).resolve().parents[2]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
import os

import numpy as np
import pandas as pd
from time import time
import shutil
import multiprocessing as mp
# Add project root to path
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))
# Import your modules directly from src_python
sys.path.insert(0, str(project_root / "src" / "src_python"))
from src.src_python.Direct_model import setup_ginette
from src.src_python.Init_folders import compile_ginette
# Get current script directory
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = str(SCRIPT_DIR)
# =============================================================================
# MODEL GEOMETRY AND DISCRETIZATION SETUP
# =============================================================================
# Set up directories and data paths
# This path contains the observational data (temperature sensors, pressure measurements)
Obs_data = os.path.join(BASE_DIR,'OBS_point/Point3_540_SOULTZ/')
# Start date and time for the simulation
# This corresponds to when field measurements began
date_simul_bg = pd.to_datetime("2022/04/21 14:00:00")


# TEMPORAL DISCRETIZATION
# Time step in seconds (900s = 15 minutes)
# This matches the measurement frequency and ensures numerical stability
dt = 900
nb_day = 30      # Simulation duration in days
# Simulation state configuration:
# 0 = steady state (time-independent, equilibrium conditions)
# 1 = transient state (time-dependent, dynamic evolution)
# We use transient state to capture temporal variations in temperature and flow
state = 1
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
# HYDROGEOLOGICAL PARAMETER DEFINITION
# =============================================================================

# GEOLOGICAL HETEROGENEITY
# Number of geological facies (material zones) in the column
# nb_zone = 1: Homogeneous porous medium
# nb_zone = 2: Two-layer system (typical for streambed environments)
nb_zone = 1

# Boundary between geological zones [m]
# Negative value indicates depth below streambed surface
alt_thk = -0.32  # Interface at 32 cm depth

# Initialize Ginette model files and return observation depths
# This function creates all necessary input files for the Ginette model:
# - Parameter files (E_parametre.dat)
# - Thermal parameter files (E_p_therm.dat)

print(f"Model domain setup:")
print(f"- Vertical extent: {z_top} to {z_bottom} m")
print(f"- Cell size: {dz} m")
print(f"- Number of cells: {int(az/dz)}")
print(f"- Time step: {dt} s ({dt/60} minutes)")
print(f"- Observation depths: {dz_obs} m")



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
    from Direct_model import (setup_ginette2,
                               initial_conditions,
                               boundary_conditions,
                               run_direct_model,
                               smooth_square_wave)
except Exception:
     # fallback to legacy module name if present
    from direct_model_ginette import (setup_ginette2,
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

def run_ginette(ID, k, n,lam,c,z_top, z_bottom, dz, dt, state, nb_day, dz_obs, date_simul_bg,nb_zone,alt_thk):
    # Temp dir:
    temp_dir = os.path.join(BASE_APP_DIR, "temp", f"temp_{ID}")
    os.makedirs(temp_dir, exist_ok=True)
    prepare_ginette_directories(temp_dir)

    def _copy_from_app(name):
        candidates = [
            os.path.join(BASE_APP_DIR, "GINETTE_SENSI", name),
            os.path.join(BASE_APP_DIR, name),
            os.path.join(REPO_ROOT, "application", "1D_Stream_aquifer_GridSearch", "GINETTE_SENSI", name),
            os.path.join(REPO_ROOT, "application", "1D_Stream_aquifer_GridSearch", name),
        ]
        for src in candidates:
            if os.path.exists(src):
                copy_file(src, temp_dir)
                return
        raise FileNotFoundError(f"Template '{name}' not found. Looked in: {candidates}")

    # Copy required template files from the application folder (tries multiple locations)
    for fname in [
        "E_parametre_bck.dat",
        "E_p_therm_bck.dat",
        "E_cdt_aux_limites_bck.dat",
        "E_zone_parameter_bck.dat",
        "ginette",
        "Sim_temperature_maille1_t.dat",
        "Sim_temperature_maille2_t.dat",
        "Sim_temperature_maille3_t.dat",
        "E_cdt_initiale.dat",
        "E_charge_initiale.dat",
        "E_cdt_aux_limites.dat",
        "E_charge_t.dat",
        "E_temp_t.dat",
        "E_colonne.dat", 
        "E_coordonnee.dat",
          "E_def_maille.dat",
          "E_temperature_initiale.dat",
            "E_zone.dat"
    ]:
        _copy_from_app(fname)

    # rename backup parameter file if present
    src_param = os.path.join(temp_dir, "E_parametre_bck.dat")
    if os.path.exists(src_param):
        os.rename(src_param, os.path.join(temp_dir, "E_parametre_backup.dat"))

    # run inside the temp directory (use absolute path)
    os.chdir(temp_dir)
    print("Current working directory:", os.getcwd())


    os.chdir(os.path.join(os.getcwd(), temp_dir))
    print(os.getcwd())
    z_obs = setup_ginette2(dt, state, nb_day, z_top, z_bottom, az, dz,
                           date_simul_bg, dz_obs)


    


    # 4) RUN SIMULATION
    sim_temp = run_direct_model(date_simul_bg,
                                z_bottom,
                                dz,
                                nb_zone,
                                alt_thk,
                                k,
                                REF_n,
                                lam,
                                REF_r,
                                REF_k2=None,
                                REF_n2=None,
                                REF_l2=None,
                                REF_r2=None)

    # Save results:
    os.chdir(os.path.join("..", ".."))

    sim_temp.to_csv(os.path.join("results",
                                 f"sim_temp_{ID}.txt"), sep=" ")

    # Del temp
    shutil.rmtree(temp_dir)
    return




if __name__ == "__main__":

    grid = pd.read_csv(os.path.join("grid_search.csv"), delimiter=";")
    os.makedirs("results", exist_ok=True)
    os.makedirs("temp", exist_ok=True)
    compile_ginette()
    
    # - Grid coordinates and observation points
    z_obs = setup_ginette(dt, state, nb_day, z_top, z_bottom, az, dz, date_simul_bg, dz_obs)

    
#    done = sorted([int(f.split("_")[-1].split(".")[0])
#                   for f in os.listdir(os.path.join("results"))])[1:]
#    remains = [i for i in grid.ID if i not in done]

#    params = [[r.ID, np.log10(r.log_k), r.lam] for r in grid.itertuples()
#              if r.ID in remains]
#    to = time()
#    with mp.Pool(processes=mp.cpu_count()-2) as pool:
#        result = pool.starmap_async(run_ginette, params)
#        pool.close()
#        pool.join()
 #   tf = time()
 #   print(f"Run time for {len(params)} simulations: {round(tf-to, 2)} s")
