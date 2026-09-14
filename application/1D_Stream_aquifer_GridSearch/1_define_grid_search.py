#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 20:19:17 2024

@author: Agnès Rivière, Samuel Larance

Construit la table de combinaisons de paramètres (grid_search.csv) utilisée
pour le calage de Ginette. Le nombre de zones géologiques (NB_ZONE) et la
liste des paramètres réellement testés (Name_parameters) sont les deux choses
à changer pour adapter ce script à un nouveau cas : voir PARAMÈTRES UTILISATEUR
ci-dessous.
"""

# %% IMPORTS
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[2]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
import os
import itertools
import numpy as np
import pandas as pd

sys.path.insert(0, str(project_root / "src" / "src_python"))

# %% =============================================================================
# PARAMÈTRES UTILISATEUR
# =============================================================================
# NB_ZONE, ALT_THK, N, Name_parameters, les bornes [min, max] de chaque
# paramètre et PARAM_STEP vivent tous dans config_lomos.py (un seul endroit à
# changer, partagé avec 2_run_real_case.py/3_misfit.py/4_plot_results.py).
sys.path.insert(0, str(project_root / "application" / "1D_Stream_aquifer_GridSearch"))
from config_lomos import NB_ZONE, ALT_THK, N, Name_parameters, PARAM_STEP, log_k, lam, n, c

NB_parameters = len(Name_parameters)
if NB_ZONE == 2 and not any(p.endswith("2") for p in Name_parameters):
    print("Avertissement : NB_ZONE=2 mais aucun paramètre de zone 2 "
          "(ex: log_k2, lam2, n2) n'est dans Name_parameters — la zone 2 "
          "utilisera des valeurs de référence fixes, pas calibrées.")

# Vérifie que chaque nom de Name_parameters correspond bien à une variable définie
for param in Name_parameters:
    if param not in globals():
        raise ValueError(f"Parameter {param} is not defined in the script.")
# =============================================================================


# Find project-relative src and application directories (no absolute paths)
def find_project_paths(start_file=__file__):
    p = Path(start_file).resolve().parent
    repo_root = None
    src_py = None
    app_dir = None
    for _ in range(8):
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

if SRC_PY.is_dir() and str(SRC_PY) not in sys.path:
    sys.path.insert(0, str(SRC_PY))

# Un sous-dossier par point (results/{POINT_NAME}/, config_lomos.py partagé
# avec les autres scripts) : un autre point lancé ensuite n'écrase pas la
# grille/les résultats de celui-ci.
sys.path.insert(0, str(BASE_APP_DIR))
from config_lomos import POINT_NAME
RESULTS_DIR = BASE_APP_DIR / "results" / POINT_NAME
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

if SRC_PY.is_dir():
    try:
        print(f"Using src_python: {SRC_PY}")
        print("Files:", sorted(os.listdir(SRC_PY)))
    except Exception:
        pass

# %% CREATE TABLE FOR THE GRID SEARCH
# %% CREATE GRIDS FOR EACH PARAMETER (générique : lit les bornes via globals(),
# donc fonctionne sans modification pour 1, 2, 3 ou 4 paramètres dans
# Name_parameters - il suffit d'éditer cette liste, rien d'autre à changer.
# Pour chaque paramètre :
#   - listé dans PARAM_STEP -> grille à pas fixe (np.arange + round)
#   - sinon [min, max] (2 valeurs) -> N valeurs régulières (np.linspace)
#   - sinon (1 valeur fixe, ou liste explicite de valeurs à tester) -> valeurs
#     utilisées telles quelles, sans balayage (comme c_grid=c chez l'étudiant)
param_grids = {}
for param_name in Name_parameters:
    param_range = globals()[param_name]
    if param_name in PARAM_STEP:
        step, decimals = PARAM_STEP[param_name]
        param_grids[param_name] = np.round(
            np.arange(param_range[0], param_range[1] + step, step), decimals)
    elif len(param_range) == 2:
        param_grids[param_name] = np.linspace(param_range[0], param_range[1], N)
    else:
        param_grids[param_name] = np.array(param_range)
    print(f"Grid for {param_name}: {param_grids[param_name]}")

# %% CREATE GRID SEARCH TABLE DYNAMICALLY
grid_search_list = list(itertools.product(*[param_grids[p] for p in Name_parameters]))
grid_search = pd.DataFrame(grid_search_list, columns=Name_parameters)

# %% remove rows with identical parameters (if any)
grid_search = grid_search.drop_duplicates().reset_index(drop=True)

# %% ADD ID, NB_ZONE AND ALT_THK COLUMNS (pour que 2_run_simulations.py sache
# quelle configuration de zones a servi à générer cette table)
grid_search.insert(0, "ID", range(len(grid_search)))
grid_search["nb_zone"] = NB_ZONE
grid_search["alt_thk"] = ALT_THK

# %% REMOVE DUPLICATES (IF ANY)
grid_search = grid_search.drop_duplicates().reset_index(drop=True)

# %% SAVE GRID SEARCH TABLE
grid_search.to_csv(RESULTS_DIR / "grid_search.csv", sep=";", index=False)
print(f"\nGrid search table saved with {len(grid_search)} combinations")
print(grid_search.head())
