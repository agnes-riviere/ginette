#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 20:19:17 2024

@author: Maxime GAUTIER,Agnes Rivière, Samuel Larance
"""


# %% IMPORTS:

import sys
from pathlib import Path
# Add project root to path
project_root = Path(__file__).resolve().parents[2]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
import os
import itertools
import importlib
import os
import numpy as np
import pandas as pd
from time import time
import shutil
import multiprocessing as mp

import glob
# Import your modules directly from src_python
sys.path.insert(0, str(project_root / "src" / "src_python"))


# %% PARAMETERS:

# GEOLOGICAL HETEROGENEITY
# Number of geological facies (material zones) in the column
# nb_zone = 1: Homogeneous porous medium
# nb_zone = 2: Two-layer system (typical for streambed environments)

# Boundary between geological zones [m]
# Negative value indicates depth below streambed surface
#alt_thk = -0.32  # Interface at 32 cm depth

# =============================================================================
# ZONE 1 PARAMETERS (Upper layer: 0 to -32 cm)
# Typically represents looser, more permeable streambed sediments
# =============================================================================

# POROSITY (φ) - Fraction of void space available for fluid flow
# High porosity typical of loose gravel/sand

# INTRINSIC PERMEABILITY (k) - Measure of medium's ability to transmit fluid
# Relationship: k = K·μ/(ρ·g) where:
# - K: hydraulic conductivity [m/s]
# - μ: dynamic viscosity [Pa·s]
# - ρ: fluid density [kg/m³]
# - g: gravitational acceleration [m/s²]
# Log₁₀(permeability in m²)
               # k = 10^(-13.5) ≈ 3.16×10^(-14) m²
               # Corresponds to coarse sand/fine gravel

# THERMAL CONDUCTIVITY (λ) - Heat conduction efficiency [W/m·K]
2.0    # Typical for water-saturated sediments

# SOLID GRAIN DENSITY (ρₛ) - Density of mineral grains [kg/m³]
# Typical for quartz-rich sediments

# HEAT CAPACITY CALCULATION (handled internally by Ginette):
# c_pm = c_w·ρ_w·n·S + c_s·ρ_s·(1-n) + c_a·ρ_a·n·(1-S)
# where:
# - c_w = 4185 J/kg·°C (specific heat of water)
# - c_s = variable J/kg·°C (specific heat of solid)
# - S = saturation ratio
# - Fixed ρ_m = 1000 kg/m³ (imposed by Ginette)

# =============================================================================
# ZONE 2 PARAMETERS (Lower layer: -32 to -40 cm)
# Typically represents more compact, less permeable substrate
# =============================================================================

#if nb_zone == 2:
#    REF_n2 = 0.3    # Lower porosity (more compact sediment)
#    REF_k2 = -13     # Higher permeability: k = 10^(-13) m²
#    REF_l2 = 2.65    # Higher thermal conductivity (more mineral content)
#    REF_r2 = 3500    # Same grain density (similar mineral composition)
# Define ranges for test parameters in the grid search:
N = 4  # number of test by parameter (nb)
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
