#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 14:00:41 2025

@author: Maxime GautIER, Agnes  Riviere, Samuel Larance

"""


# %% IMPORTS:
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import sys
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
                               smooth_square_wave,remove_first_two_days_time_based)
except Exception:
    print("Error importing Direct_model from src_python.")

try:
    from Init_folders import prepare_ginette_directories
except Exception:
    print("Error importing Init_folders from src_python.")

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



# %% FUNCTIONS:
def misfit(obs, sim, err):
    return np.sum(((sim - obs) / err)**2)


def likelyhood(m):
    return np.exp(-0.5*m)


def log_likelyhood(m):
    return np.log(likelyhood(m))

# %% MISFIT:
# in repertory results, for each simulation done, compute misfit with observed data
results = pd.read_csv(os.path.join(RESULTS_DIR,"grid_search.csv"), delimiter=";")
obs_data = pd.read_csv(os.path.join(RESULTS_DIR,"observed_data.txt"), delimiter=" ",
                       index_col=[0])



results["misfit_1"] = np.nan
results["misfit_2"] = np.nan
results["misfit_3"] = np.nan

done = sorted([
    int(f.split("_")[-1].split(".")[0])
    for f in os.listdir(RESULTS_DIR)
    if f.startswith("sim_temp_") and f.endswith('.txt')
])

print("dans done il y a :", done)

# Misfit by simulations:
err = 0.2
results["err_misfit"] = err
for i in tqdm(done, desc="Compute misfit"):
    sim = pd.read_csv(os.path.join(RESULTS_DIR, f"sim_temp_{i}.txt"),
                      delimiter=" ", index_col=[0])
    remove_first_two_days_time_based(sim, obs_data)
    m1 = misfit(obs_data.Temp1, sim.Temp1, err=err)
    m2 = misfit(obs_data.Temp2, sim.Temp2, err=err)
    m3 = misfit(obs_data.Temp3, sim.Temp3, err=err)

    results.loc[results["ID"] == i, "misfit_1"] = m1
    results.loc[results["ID"] == i, "misfit_2"] = m2
    results.loc[results["ID"] == i, "misfit_3"] = m3

# Total misfit and likelihood
results["misfit_tot"] = (results["misfit_1"] + results["misfit_2"] +
                         results["misfit_3"])

results["L1"] = likelyhood(results["misfit_1"])
results["L2"] = likelyhood(results["misfit_2"])
results["L3"] = likelyhood(results["misfit_3"])

results.to_csv(os.path.join(RESULTS_DIR,"results.txt"), sep=" ")