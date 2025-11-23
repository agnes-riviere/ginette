#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 13:43:32 2025

@author: mgautier
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 09:51:06 2025

@author: Maxime Gautier
@modif: Agnès Rivière, Samuel Larance
"""


# IMPORT:
import os
import sys
import numpy as np
import pandas as pd
from time import time
import shutil
import multiprocessing as mp
from pathlib import Path

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


def run_ginette(ID, k, lam):
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
    ]:
        _copy_from_app(fname)

    # rename backup parameter file if present
    src_param = os.path.join(temp_dir, "E_parametre_bck.dat")
    if os.path.exists(src_param):
        os.rename(src_param, os.path.join(temp_dir, "E_parametre_backup.dat"))

    # run inside the temp directory (use absolute path)
    os.chdir(temp_dir)
    print("Current working directory:", os.getcwd())

    # Domain paramters:
    z_top = 0
    z_bottom = -5
    n_depths = 250
    dz_obs = 0.1
    az = abs(z_top - z_bottom)
    dz = az / n_depths

    # Time parameters:
    nb_day = 30
    dt = 600
    state = 1  # 1 for transcient
    start_date = "2022/04/21 14:00:00"
    date_simul_bg = pd.to_datetime(start_date)

    print("Setup ginette model...")
    z_obs = setup_ginette2(dt, state, nb_day, z_top, z_bottom, az, dz,
                           date_simul_bg, dz_obs)
    # Moded parameters:
    nb_zone = 1
    alt_thk = 0
    REF_n = 0.05
    REF_r = 3500

    # 1) PREPARE BOUNDARY CONDITIONS:
    # Define data table:
    seconds_per_day = 86400
    seconds_per_week = 604800
    t_final = nb_day * seconds_per_day
    df_BC = pd.DataFrame()
    df_BC["times"] = np.arange(0, t_final, dt)
    df_BC["days"] = df_BC["times"]/seconds_per_day

    # Calculate boundaries temperature:
    deg_per_day = 0.3
    df_BC["T_top"] = (5*np.sin(2*np.pi*df_BC["times"]/seconds_per_day) +
                      deg_per_day*df_BC["days"] +
                      3*np.sin(2*np.pi*df_BC["times"] / seconds_per_week) + 17)
    df_BC["T_bottom"] = 10 * np.ones_like(df_BC["times"])

    # Calculate hydraulic head at top and bottom boundaries:
    reduced_period = 3 * seconds_per_week
    base_head = 5.0
    head_amplitude = 0.2
    df_BC["h_top"] = (head_amplitude * smooth_square_wave(
        df_BC["times"], reduced_period) + base_head)
    df_BC["h_bottom"] = base_head * np.ones_like(df_BC["times"])
    df_BC["head_gradient"] = (df_BC["h_top"] - df_BC["h_bottom"]) / az

    # 2) CREATE MEASUREMENT TABLE:
    obs_temp = pd.DataFrame({
        "dates": date_simul_bg + pd.to_timedelta(df_BC["times"], unit="s"),
        "h_top": df_BC["h_top"],
        "h_bottom": df_BC["h_bottom"],
        "T_top": df_BC["T_top"],
        "T_bottom": df_BC["T_bottom"]})

    # 3) SET INITIAL ANd BOUNDARIES CONDITIONS:
    # Set initial conditions:
    z_obs = [-5]
    initial_conditions(obs_temp, z_top, z_bottom, dz, z_obs)

    # Set boundary conditions in temperature:
    boundary_conditions(obs_temp, dt)

    # Save initial and boundary conditions file in /home/ariviere/Programmes/ginette/application/1D_Stream_aquifer_GridSearch/SYNTHETIC_CASES:
    # E_charge_initiale.dat, E_charge_t.dat, E_temp_t.dat, E_temperature_initiale.dat
    for fname in [ 
        "E_charge_t.dat",
        "E_temp_t.dat"
    ]:
        src = os.path.join(temp_dir, fname)
        dst_dir = os.path.join(BASE_APP_DIR, "SYNTHETIC_CASES")
        if os.path.exists(src):
            copy_file(src, dst_dir)
        else:
            print(f"Warning: {fname} not found in {temp_dir}")
        

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
    sim_temp.to_csv(os.path.join(RESULTS_DIR,
                                 f"sim_temp_{ID}.txt"), sep=" ")

    # Del temp
    shutil.rmtree(temp_dir)
    return

run_ginette(-1, -12, 2.5)

# %% PLOT OBSERVED DATA:
data = pd.read_csv(os.path.join(RESULTS_DIR, "sim_temp_-1.txt"), delimiter=" ",
                   index_col=[0])
data.to_csv(os.path.join(RESULTS_DIR, "observed_data.txt"), sep=" ")




import matplotlib.pyplot as plt
plt.rcParams["font.size"] = 16

fig, ax = plt.subplots(figsize=(16, 9), dpi=100)
ax.plot(data.Time, data.Temp1, ls="-", lw=2, color="red",
        label="Temperature 1 (-10 cm)")
ax.plot(data.Time, data.Temp2, ls="-", lw=2, color="green",
        label="Temperature 2 (-20 cm")
ax.plot(data.Time, data.Temp3, ls="-", lw=2, color="blue",
        label="Temperature 3 (-30 cm)")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Temperature (°C)")
ax.grid()
ax.legend(loc="upper right")
ax.set_title("Observed data", loc="left")

# %%
