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
n = 0.05
c= 3500
k=-12.5
lam=2.5

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
# For 4 parameters, we can plot 2 marginals 1D and 1 2D plot
# plus the temperature time series for the best model   

# %% PLOT: Adaptive for any number of parameters (NB_parameters, Name_parameters)


# Determine grid layout based on number of parameters
if NB_parameters <= 2:
    grid_rows, grid_cols = 2, 2
elif NB_parameters <= 4:
    grid_rows, grid_cols = 2, 2
else:
    grid_rows, grid_cols = 3, 2

fig = plt.figure(figsize=(16, 9 * (grid_rows / 2)), dpi=100)
gs = gridspec.GridSpec(grid_rows, grid_cols, figure=fig)

# Plot marginals (one per subplot, left column)
axes_list = []
for idx, param in enumerate(Name_parameters):
    if idx < grid_rows:
        ax = fig.add_subplot(gs[idx, 0])
        ax.plot(marginals[param].index, marginals[param].values, lw=2, ls="-",
                color="black", label="Marginal 1D")
        ax.axvline(results.loc[best, param], lw=3, ls="--", color="gold",
                   label="Best value")
        ax.set_xlabel(param)
        ax.set_ylabel(assess_var)
        ax.grid()
        ax.legend(loc="upper left", fontsize=10)
        axes_list.append(ax)

# Plot 2D scatter (first two parameters, right column, top)
param1, param2 = Name_parameters[0], Name_parameters[1]
axb = fig.add_subplot(gs[0, 1])
gci = axb.scatter(results[param1], results[param2], c=results.misfit_tot, cmap=cmap, s=50)
fig.colorbar(gci, ax=axb, label=assess_var)
axb.set_xlabel(param1)
axb.set_ylabel(param2)
axb.grid()
axb.scatter(results.loc[best, param1], results.loc[best, param2],
            marker="*", edgecolors="black", facecolors="gold", s=200, lw=1, label="Best")
axb.legend(loc="upper left")

# Plot temperatures (right column, bottom)
axc = fig.add_subplot(gs[1, 1])
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
axc.legend(loc="lower right", fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "grid_search_results.png"), dpi=150, bbox_inches="tight")
plt.show()




# %% PLOT: Create a figure for each pair of parameters

from itertools import combinations

# Get all pairs of parameters
param_pairs = list(combinations(Name_parameters, 2))

for pair_idx, (param1, param2) in enumerate(param_pairs):
    # Determine grid layout based on number of parameters
    if NB_parameters <= 2:
        grid_rows, grid_cols = 2, 2
    elif NB_parameters <= 4:
        grid_rows, grid_cols = 2, 2
    else:
        grid_rows, grid_cols = 3, 2

    fig = plt.figure(figsize=(16, 9 * (grid_rows / 2)), dpi=100)
    gs = gridspec.GridSpec(grid_rows, grid_cols, figure=fig)

    # Plot marginals for all parameters (left column)
    axes_list = []
    for idx, param in enumerate(Name_parameters):
        if idx < grid_rows:
            ax = fig.add_subplot(gs[idx, 0])
            ax.plot(marginals[param].index, marginals[param].values, lw=2, ls="-",
                    color="black", label="Marginal 1D")
            ax.axvline(results.loc[best, param], lw=3, ls="--", color="gold",
                       label="Best value")
            ax.set_xlabel(param)
            ax.set_ylabel(assess_var)
            ax.grid()
            ax.legend(loc="upper left", fontsize=10)
            axes_list.append(ax)

    # Plot 2D scatter for current pair (right column, top)
    axb = fig.add_subplot(gs[0, 1])
    gci = axb.scatter(results[param1], results[param2], c=results.misfit_tot, cmap=cmap, s=50)
    fig.colorbar(gci, ax=axb, label=assess_var)
    axb.set_xlabel(param1)
    axb.set_ylabel(param2)
    axb.grid()
    axb.scatter(results.loc[best, param1], results.loc[best, param2],
                marker="*", edgecolors="black", facecolors="gold", s=200, lw=1, label="Best")
    axb.legend(loc="upper left")
    axb.set_title(f"2D: {param1} vs {param2}")

    # Plot temperatures (right column, bottom)
    axc = fig.add_subplot(gs[1, 1])
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
    axc.legend(loc="lower right", fontsize=9)

    plt.tight_layout()
    # Save with pair name in filename
    filename = f"grid_search_results_{param1}_vs_{param2}.png"
    plt.savefig(os.path.join(RESULTS_DIR, filename), dpi=150, bbox_inches="tight")
    print(f"Saved: {filename}")
    plt.show()