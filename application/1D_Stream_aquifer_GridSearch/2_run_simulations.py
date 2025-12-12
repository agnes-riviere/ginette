#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 09:51:06 2025
For the SYNTHETIC_CASES defined in the application folder,
@author: Maxime Gautier, Agnes Rivière, Samuel Larance
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
# Import your modules directly from src_python
sys.path.insert(0, str(project_root / "src" / "src_python"))



Delete_sim="True"


def delete_sim_temp(my_dir,temp):
    """
    Delete files in my_dir whose name starts with 'sim_temp'.
    Deletes temp_ repertory in temp if exists.
    -----------
    my_dir : str
        Directory path where to delete files.
    -----------
    Returns
    None
    ----------- 
    """
    p = Path(my_dir)
    if not p.exists():
        return
    if not p.is_dir():
        raise NotADirectoryError(f"Not a directory: {my_dir}")
    for f in p.glob("sim_temp*"):
        try:
            if f.is_file():
                f.unlink()
        except Exception as e:
            print(f"Failed to remove {f}: {e}")
    # delete temp_ repertory in temp if exists
    temp_p = Path(temp)
    if temp_p.exists() and temp_p.is_dir():
        for sub in temp_p.glob("temp_*"):
            try:
                if sub.is_dir():
                    shutil.rmtree(sub)
            except Exception as e:
                print(f"Failed to remove directory {sub}: {e}")



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
                               smooth_square_wave,generate_zone_parameters)
except Exception:
     # fallback to legacy module name if present
    from direct_model_ginette import (setup_ginette2,
                                       initial_conditions,
                                       boundary_conditions,
                                       run_direct_model,
                                       smooth_square_wave)


try:
    from Init_folders import prepare_ginette_directories,compile_ginette_src
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

def run_ginette(ID, k, n,lam,c):
    # Temp dir:
    temp_dir = os.path.join(BASE_APP_DIR, "temp", f"temp_{ID}")
    os.makedirs(temp_dir, exist_ok=True)


    def _copy_from_app(name):
        candidates = [
            os.path.join(BASE_APP_DIR, "SYNTHETIC_CASES", name),
            os.path.join(BASE_APP_DIR, name),
            os.path.join(REPO_ROOT, "application", "1D_Stream_aquifer_GridSearch", "SYNTHETIC_CASES", name),
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
        "E_cdt_initiale.dat",
    ]:
        _copy_from_app(fname)



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
    nb_zone=1
    alt_thk=0
    # Time parameters:
    nb_day = 30
    dt = 600
    state = 1  # 1 for transcient
    start_date = "2022/04/21 14:00:00"
    date_simul_bg = pd.to_datetime(start_date)

    print("Setup ginette model...")
    z_obs = setup_ginette2(dt, state, nb_day, z_top, z_bottom, az, dz,
                           date_simul_bg, dz_obs)


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
    #compile ginette fortran source code from GINETTE src folder
    compile_ginette_src(REPO_ROOT)
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
                                n,
                                lam,
                                c,
                                REF_k2=None,
                                REF_n2=None,
                                REF_l2=None,
                                REF_r2=None)



    # Save results:
    os.chdir(os.path.join("..", ".."))
    print('on est ici',os.getcwd())
    sim_temp.to_csv(os.path.join(RESULTS_DIR,
                                 f"sim_temp_{ID}.txt"), sep=" ")

    # Del temp
     #    shutil.rmtree(temp_dir)
    return


if __name__ == "__main__":

    # result is ginette/application/1D_Stream_aquifer_GridSearch/results/
    # temp is ginette/application/1D_Stream_aquifer_GridSearch/temp/
    # verify repertory temp and results exist
    os.makedirs(os.path.join(BASE_APP_DIR, "temp"), exist_ok=True)
    os.makedirs(os.path.join(BASE_APP_DIR, "results"), exist_ok=True)
    if (Delete_sim=="True"):
        delete_sim_temp(os.path.join(BASE_APP_DIR, "results"),os.path.join(BASE_APP_DIR, "temp"))
    # grid search in results/grid_search.csv
    grid = pd.read_csv(os.path.join(BASE_APP_DIR, "results", "grid_search.csv"), delimiter=";")

    # find which simulations are already done
    # if file sim_temp_ID.txt exists in results/, consider it done
    if (Delete_sim!="True"):
        done = sorted([int(f.split("_")[-1].split(".")[0])
                   for f in os.listdir(os.path.join(BASE_APP_DIR, "results"))])[1:]
        remains = [i for i in grid.ID if i not in done]
    else:
        remains = grid.ID.tolist()

    params = [[r.ID, r.log_k,r.poro, r.lam,r.cap] for r in grid.itertuples()
              if r.ID in remains]
    to = time()
    with mp.Pool(processes=mp.cpu_count()-2) as pool:
        result = pool.starmap_async(run_ginette, params)
        pool.close()
        pool.join()
    tf = time()
    print(f"Run time for {len(params)} simulations: {round(tf-to, 2)} s")
