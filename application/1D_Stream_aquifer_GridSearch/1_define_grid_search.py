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
# NOMBRE DE ZONES GÉOLOGIQUES DE LA COLONNE
#   1 = milieu homogène
#   2 = deux couches séparées par une interface à ALT_THK. NB : dans
#       generate_zone_parameters, la ZONE 2 est la couche DU HAUT (entre la
#       surface et ALT_THK) et la ZONE 1 est la couche DU BAS (entre ALT_THK et
#       le fond) — voir coord['zone']=np.where(z>=alt_thk,2,1) dans Direct_model.py.
# NB_ZONE=1 retenu suite au diagnostic TempMolo (sur-amplification de Temp1,
# non liée à une hétérogénéité de matériau) - historique complet dans la
# mémoire projet "TempMolo mal positionné lomos230/231".
NB_ZONE = 1

# Profondeur de l'interface [m], négative = sous la surface du lit du ru.
# Sans effet tant que NB_ZONE=1 (conservée pour repasser à 2 zones facilement).
ALT_THK = -0.05

# Nombre de valeurs testées par paramètre : sert uniquement de repli
# (np.linspace) pour tout paramètre absent de PARAM_STEP ci-dessous.
N = 8

# Pas fixe (résolution) par paramètre, comme dans 1_define_grid_search_tris.py
# (script de l'étudiant) : {nom: (pas, décimales d'arrondi)}. Remplace le N
# unique ci-dessus pour les paramètres listés ici - la grille est alors
# np.arange(min, max+pas, pas) (le +pas inclut la borne max, que arange exclut
# sinon), arrondie à `décimales` pour éviter les valeurs parasites du type
# 0.300000000000004 dues aux erreurs de floating point de np.arange. Un
# paramètre de Name_parameters absent d'ici retombe sur N/np.linspace.
PARAM_STEP = {
    "log_k": (1, 1),
    "lam": (1, 1),
    "n": (0.05, 3),
}

# Parameters ranges [min, max]:
log_k = [-18,-11]  #  # log10(perméabilité intrinsèque k [m2])
lam = [1,6]   # conductivité thermique de la fraction SOLIDE [W/m/K]
n=[5e-2, 0.6 ]  # porosité [-]
c=[2000]   # capacité thermique volumique de la fraction solide [J/m3/K]
# PARAMÈTRES RÉELLEMENT TESTÉS DANS LE GRID SEARCH : chaque nom doit correspondre
# à une variable [min, max] définie ci-dessus. Ajouter ou retirer un paramètre ne
# demande de modifier que cette liste (la construction de la grille plus bas est
# générique et lit les bornes via globals()).
Name_parameters = ["log_k", "lam", "n", "c"]

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

RESULTS_DIR = BASE_APP_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

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
