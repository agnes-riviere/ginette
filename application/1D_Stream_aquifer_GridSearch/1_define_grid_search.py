#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 20:19:17 2024

@author: Maxime GAUTIER,Agnes Rivière, Samuel Larance
"""


# %% IMPORTS:
import os
import itertools
import numpy as np
import pandas as pd
import sys
from pathlib import Path


# %% PARAMETERS:
# Define ranges for test parameters in the grid search:
N = 3  # number of test by parameter (nb)
NB_parameters = 4  # number of parameters to test (nb)
Name_parameters = ["log_k", "lam", "poro","cap"]  # names of parameters to test (list of str)
# check NB_parameters and Name_parameters consistency
if NB_parameters != len(Name_parameters):
    raise ValueError("NB_parameters and Name_parameters length mismatch.")

"""
# Example values:
log_k = -12  # log of intrinsec permeability (k in m^2)
lam = 2.5  # thermal conductivity (W/m.°C)
"""

# Parameters ranges [min, max]:
log_k = [-15, -12]  # log of intrinsec permeability (k in m^2)
lam = [1, 3]  # thermal conductivity (W/m.°C)
poro=[0.1, 0.5]  # porosity (-)
cap=[1.2e3, 3.2e3]  # volumetric heat capacity (J/m^3.°C)

# Check all parameter intervals are all defined and exist 
# i.e all name_parameters are defined in a vector for example log_k, lam, poro, cap
# je ne veux pas que cela soit en dur les noms peuvent changer
for param in Name_parameters:
    if param not in globals():
        raise ValueError(f"Parameter {param} is not defined in the script.")
 

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




# %% CREATE TABLE FOR THE GRID SEARCH
# %% more generaliste take the name of parameters and their ranges from a dictionnary
# %% CREATE GRIDS FOR EACH PARAMETER
param_grids = {}
for param_name in Name_parameters:
    param_range = globals()[param_name]
    param_grids[param_name] = np.linspace(param_range[0], param_range[1], N)
    print(f"Grid for {param_name}: {param_grids[param_name]}")

# %% CREATE GRID SEARCH TABLE DYNAMICALLY
grid_search_list = list(itertools.product(*[param_grids[p] for p in Name_parameters]))
grid_search = pd.DataFrame(grid_search_list, columns=Name_parameters)

# %% remove rows with identical parameters (if any)
grid_search = grid_search.drop_duplicates().reset_index(drop=True)

# %% ADD ID COLUMN
grid_search.insert(0, "ID", range(len(grid_search)))

# %% REMOVE DUPLICATES (IF ANY)
grid_search = grid_search.drop_duplicates().reset_index(drop=True)

# %% SAVE GRID SEARCH TABLE
grid_search.to_csv(os.path.join(RESULTS_DIR, "grid_search.csv"), sep=";", index=False)
print(f"\nGrid search table saved with {len(grid_search)} combinations")
print(grid_search.head())
