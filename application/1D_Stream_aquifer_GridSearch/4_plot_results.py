#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 15:21:53 2025

@author: Maxime GAUTIER, Agnes Rivière, Samuel Larance
"""


# %% IMPORTS:
import os
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.pyplot as plt
cmap = matplotlib.colormaps['nipy_spectral'].resampled(256).reversed()
fontsize = 12
plt.rcParams["font.size"] = fontsize
from tqdm import tqdm
from time import time
import shutil
import multiprocessing as mp
from pathlib import Path
import sys
# %% PARAMETERS:
# Define ranges for test parameters in the grid search:
N = 2  # number of test by parameter (nb)
NB_parameters = 4  # number of parameters to test (nb)
Name_parameters = ["log_k", "lam", "poro","cap"]  # names of parameters to test (list of str)
# check NB_parameters and Name_parameters consistency
if NB_parameters != len(Name_parameters):
    raise ValueError("NB_parameters and Name_parameters length mismatch.")
n = 0.25
c= 3500
k=-12.5
lam=2.5
REF=[k,lam,n,c]  # reference values for the parameters (list of float)
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
    from direct_model import (setup_ginette2,
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


# %% LOAD RESULTS:
results = pd.read_csv(os.path.join(RESULTS_DIR, "results.txt"), delimiter=" ",
                      index_col=[0])


assess_var = "misfit_tot"
results.sort_values(by=assess_var, inplace=True, ascending=True)
best = results[results[assess_var] == results[assess_var].min()].ID.values[0]
sim_data = pd.read_csv(os.path.join(RESULTS_DIR, f"sim_temp_{best}.txt"),
                       delimiter=" ", index_col=[0])
obs_data = pd.read_csv(os.path.join(RESULTS_DIR,"observed_data.txt"), delimiter=" ",
                       index_col=[0])

# Marginals 1D:
# Marginals 1D: marginal with all parameters in Name_parameters
marginals = {}
for param in Name_parameters:
    marginals[param] = results.groupby(param)[assess_var].sum()
    print(f"Marginal for {param}:\n{marginals[param]}\n")

# Access : marginals["log_k"], marginals["lam"], marginals["poro"], marginals["cap"])


print(f"Best parameters (ID={best}):")
for param in Name_parameters:
    print(f"  {param} = {results.loc[best, param]}")
print(f"With total misfit: {results.loc[best, assess_var]}")    

# %% PLOT: adapt for 4 parameters
# %% PLOT: Create a figure for each pair of parameters

from itertools import combinations

# Get all pairs of parameters

if NB_parameters == 2:
    grid_rows, grid_cols = 2, 2
elif NB_parameters ==3:
    grid_rows, grid_cols = 3, 3
else:
    grid_rows, grid_cols = 4,4 

fig = plt.figure(figsize=(16, 9), dpi=100)
gs = gridspec.GridSpec(grid_rows, grid_cols, figure=fig)

axa = fig.add_subplot(gs[1, 1])  # marginals 1D k
axb = fig.add_subplot(gs[0, 1]) # 2D plot k-lam
axc = fig.add_subplot(gs[1, 0]) # time series
axd = fig.add_subplot(gs[0, 0]) # marginals 1D lam
if NB_parameters==3 :
    axe = fig.add_subplot(gs[3, 3]) # marginals 1D poro
    axf = fig.add_subplot(gs[2, 1])  # 2D plot k-poro
    axg = fig.add_subplot(gs[2, 2])  # 2D plot lam-poro
elif NB_parameters==4:
    axg = fig.add_subplot(gs[3, 0])
    axh = fig.add_subplot(gs[3, 1])

# Afficher l'index d'une marginal spécifique
for param in Name_parameters:
    print(f"Index for {param}:")
    print(marginals[param].index)
    print()
axa.plot(marginals[Name_parameters[0]].index, marginals[Name_parameters[0]].values, lw=2, ls="-",color="black", label="Marginal 1D")
axa.axvline(REF[0], lw=4, ls="--", color="lime", label="True value")
axa.axvline(results.loc[best, Name_parameters[0]], lw=3, ls="--", color="gold",
            label="Best value")
axa.legend(loc="upper left")
# Name of parameter on x axis
axa.set_xlabel(Name_parameters[0])
axa.set_ylabel(r"$\Phi{}_{tot.}$")
axa.grid()

axd.plot(marginals[Name_parameters[1]].index, marginals[Name_parameters[1]].values, lw=2, ls="-",
         color="black", label="Marginal 1D")
axd.axvline(REF[1], lw=4, ls="--", color="lime", label="True value")
axd.axvline(results.loc[best, Name_parameters[1]], lw=3, ls="--", color="gold",
            label="Best value")
axd.legend(loc="upper left")
axd.set_xlabel(Name_parameters[1])
axd.set_ylabel(r"$\Phi{}_{tot.}$")
axd.grid()



# Plot 2D:
gci = axb.scatter(results.log_k, results.lam, c=results.misfit_tot,
                  cmap=cmap)
fig.colorbar(gci, ax=axb, label=r"$\Phi{}_{tot.}$")
axb.set_xlabel(r"$log_{10}(k)$")
axb.set_ylabel(r"$\lambda\;(W/m°C)$")
axb.grid()
axb.set_xlim(results.log_k.min(), results.log_k.max())
axb.set_ylim(results.lam.min(), results.lam.max())

axb.scatter(k, lam, marker="o", edgecolors="lime", facecolors="None", s=100,
            lw=5, label="True values")
axb.scatter(results.loc[best, "log_k"],
            results.loc[best, "lam"], marker="*", edgecolors="black",
            facecolors="gold", s=80, lw=1, label="Best")
axb.legend(loc="upper left")

# Plot temperatures:
axc.plot(obs_data.Time, obs_data.Temp1, ls="-", lw=4, color="red",
         label="Temperature 1 (-10 cm)")
axc.plot(obs_data.Time, obs_data.Temp2, ls="-", lw=4, color="green",
         label="Temperature 2 (-20 cm)")
axc.plot(obs_data.Time, obs_data.Temp3, ls="-", lw=4, color="blue",
         label="Temperature 3 (-30 cm)")
axc.plot(sim_data.Time, sim_data.Temp1, ls="--", lw=2, color="orange",
         label="Sim. temp. 1")
axc.plot(sim_data.Time, sim_data.Temp2, ls="--", lw=2, color="lime",
         label="Sim. temp. 2")
axc.plot(sim_data.Time, sim_data.Temp3, ls="--", lw=2, color="dodgerblue",
         label="Sim. temp. 3")
axc.set_xlabel("Time (s)")
axc.set_ylabel("Temperature (°C)")
axc.grid()
axc.legend(loc="lower right")

plt.show()






