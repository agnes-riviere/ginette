#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Today's date: 2025-12-11
@author: Agnès Riviere
Run a set of Ginette simulations for a real case of 1D stream-aquifer interaction
Repertory of simulation setup and results: application/1D_Stream_aquifer_GridSearch/GINETTE_SENSI
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
import subprocess
import multiprocessing as mp
# Add project root to path
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))
# Import your modules directly from src_python
sys.path.insert(0, str(project_root / "src" / "src_python"))
from src.src_python.Direct_model import setup_ginette,run_direct_model,setup_ginette2,generate_zone_parameters,reuse_end_in_initial,setup_ginette_perm
from src.src_python.Init_folders import compile_ginette,compile_ginette_src,prepare_ginette_directories
# Get current script directory
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = str(SCRIPT_DIR)


#==============================================================================
# SIMULATION SETTINGS
#==============================================================================
# Delete temp folder after simulation   
Delete_sim="True"
# =============================================================================
# MODEL GEOMETRY AND DISCRETIZATION SETUP
# =============================================================================
# Set up directories and data paths
# This path contains the observational data (temperature sensors, pressure measurements)
# POINT_NAME vient de config_lomos.py (partagé avec 0_boundary_conditions_real_case.py) -
# DOIT rester identique au point/à la période utilisés là-bas : les fichiers copiés
# depuis GINETTE_SENSI/ (E_temp_t.dat, E_charge_t.dat, ...) ont été générés pour CE
# point et CETTE période précis, un décalage ici les rendrait incohérents avec les
# métadonnées (date, durée) écrites dans E_parametre.dat.
sys.path.insert(0, str(SCRIPT_DIR))
from config_lomos import (POINT_NAME, RUN_GRID_SEARCH, RUN_MISFIT, RUN_PLOTS,
                           SAVE_VELOCITY_PROFILES,
                           DATE_SIMUL_BG as date_simul_bg,
                           NB_DAY as nb_day, DT as dt, MAX_WORKERS)
Obs_data = os.path.join(BASE_DIR, 'OBS_point', POINT_NAME)
# Simulation state configuration:
# 0 = steady state (time-independent, equilibrium conditions)
# 1 = transient state (time-dependent, dynamic evolution)
# We use transient state to capture temporal variations in temperature and flow
state = 1
# SPATIAL DOMAIN DEFINITION (1D vertical column)
# TOP_SENSOR/BOTTOM_SENSOR (Config A/C, voir README.md) vivent dans
# config_lomos.py, partagés avec 0_boundary_conditions_real_case.py et
# 3_misfit.py - sinon E_temp_t.dat/E_charge_t.dat ne correspondent plus à la
# géométrie régénérée ici à chaque point de la grille.
from config_lomos import Z_TOP, Z_BOTTOM, DZ_OBS
z_top, z_bottom = Z_TOP, Z_BOTTOM
az = abs(z_top - z_bottom)  # Total column height [m]

# GRID DISCRETIZATION
dz = 0.01        # Vertical cell size [m] (1 cm resolution)
                 # Fine discretization needed to capture thermal gradients

# OBSERVATION DEPTHS
dz_obs = DZ_OBS  # Spacing between temperature sensors [m], derived from SENSOR_DEPTHS






# =============================================================================
# HYDROGEOLOGICAL PARAMETER DEFINITION
# =============================================================================
# nb_zone et alt_thk ne sont PAS fixés ici : ils sont lus depuis les colonnes
# "nb_zone"/"alt_thk" de grid_search.csv (écrites par 1_define_grid_search.py),
# pour n'avoir qu'une seule source de vérité sur la configuration de zones.

# Densité de la fraction solide : n'est plus passée via E_zone_parameter.dat
# (2026-07-22 : CASE('ZHZ') dans src/ginette_V2.f90 lit désormais cpmzone à la
# place de rhomzone sur cette colonne - voir capacité calorifique ci-dessous).
# Ginette utilise donc directement rhosi de E_p_therm.dat (1180 kg/m3, fixe,
# irhomi=0 => jamais spatialisée ; sédiment riche en matière organique plutôt
# que quartz pur, baissé depuis 2650 le 2026-07-24 pour rester cohérent avec
# c=2000 J/kg/K - voir config_lomos.py), sans qu'on ait besoin de la repasser ici.

# Valeurs de repli pour log_k/lam/n/c : utilisées UNIQUEMENT si le paramètre
# correspondant n'est pas dans Name_parameters de 1_define_grid_search.py (donc
# pas de colonne dans grid_search.csv - grid search à 1, 2, 3 ou 4 paramètres).
# REF_LOG_K/REF_LAM/REF_N = meilleur point trouvé sur lomos230 (voir mémoire
# projet "TempMolo mal positionné lomos230/231", misfit_tot=2201.3) - à changer
# si on fixe un paramètre pour un autre point ou avec une autre valeur de
# référence. REF_HEAT_CAPACITY (capacité calorifique du solide [J/kg/K],
# colonne cpmzone côté Fortran pour ytest=ZHZ) = valeur par défaut de
# 1_define_grid_search.py (c=[2000]).
REF_LOG_K = -15.86
REF_LAM = 1.57
REF_N = 0.51
REF_HEAT_CAPACITY = 2000.0

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

# --- Dossier de sortie, un sous-dossier par point (results/{POINT_NAME}/) :
# un autre point lancé ensuite n'écrase pas les résultats de celui-ci. ---
RESULTS_DIR = os.path.join(BASE_APP_DIR, "results", POINT_NAME)
os.makedirs(RESULTS_DIR, exist_ok=True)

# debug listing to verify (non-blocking)
if os.path.isdir(SRC_PY):
    try:
        print(f"Using src_python: {SRC_PY}")
        print("Files:", sorted(os.listdir(SRC_PY)))
    except Exception:
        pass

# Import project modules from src_python
try:
    # preferred: modules as they appear in src/src_python
    from src.src_python.Direct_model import (setup_ginette2,
                               initial_conditions,
                               boundary_conditions,
                               run_direct_model,
                               smooth_square_wave,
                               setup_ginette_perm,
                               generate_zone_parameters,
                               reuse_end_in_initial,
                               setup_ginette)
except Exception:
                print(f"import error in Direct_model from {SRC_PY} ")


try:
    from src.src_python.Init_folders import prepare_ginette_directories
except Exception:
    print(f"import error in Init_folders from {SRC_PY} ")

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

def run_ginette(ID, k, n,lam,c,date_simul_bg,dt,nb_day,state,z_top ,z_bottom ,az ,dz ,dz_obs  ,nb_zone ,alt_thk,
                k2=None, n2=None, lam2=None, c2=None):
    # Temp dir:
    temp_dir = os.path.join(BASE_APP_DIR, "temp", f"temp_{ID}")
    os.makedirs(temp_dir, exist_ok=True)


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
        "E_cdt_aux_limites_bck.dat",
        "E_cdt_aux_limites_perm_bck.dat",
        "E_p_therm_bck.dat",
        "E_cdt_aux_limites_bck.dat",
        "E_zone_parameter_bck.dat",
        "E_zone_parameter_bck_2zones.dat",
        "ginette",
        "E_cdt_initiale.dat",
        "E_charge_initiale.dat",
        "E_cdt_aux_limites.dat",
        "E_charge_t.dat",
        "E_pression_initiale.dat",
        "E_temp_t.dat",
        "E_colonne.dat", 
        "E_coordonnee.dat",
          "E_def_maille.dat",
          "E_temperature_initiale.dat",
            "E_zone.dat"
    ]:
        _copy_from_app(fname)


    # run inside the temp directory (use absolute path)
    os.chdir(temp_dir)

    state=0
    setup_ginette_perm(dt, state, nb_day, z_top, z_bottom, az, dz, date_simul_bg, dz_obs)
    generate_zone_parameters(z_bottom, dz, nb_zone, alt_thk, k, n, lam, c,
                              REF_k2=k2, REF_n2=n2, REF_l2=lam2, REF_r2=c2)
    subprocess.call(["./ginette"])
    # reuse end of steady state in initial conditions for transient
    source_file="S_pression_charge_temperature.dat"
    destination_file="E_charge_initiale.dat"
    reuse_end_in_initial(source_file, destination_file)
    destination_file="E_temperature_initiale.dat"
    reuse_end_in_initial(source_file, destination_file)
    state=1
    z_obs=setup_ginette(dt, state, nb_day, z_top, z_bottom, az, dz,
                 date_simul_bg, dz_obs)
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
                                REF_k2=k2,
                                REF_n2=n2,
                                REF_l2=lam2,
                                REF_r2=c2)

    # Vitesse de Darcy réellement calculée par Ginette (résolution numérique
    # complète, pas juste k*gradient) - sauvée pour le diagnostic de Péclet de
    # 3_misfit.py (qui ne s'en sert que pour la meilleure simulation, avec repli
    # analytique sinon) plutôt que de la recalculer après coup. SAVE_VELOCITY_PROFILES
    # (config_lomos.py) : False par défaut, ~9.5 Mo/simulation sinon.
    if SAVE_VELOCITY_PROFILES:
        vel_src = os.path.join(os.getcwd(), "Sim_velocity_profil_t.dat")
        if os.path.exists(vel_src):
            shutil.copy(vel_src, os.path.join(RESULTS_DIR, f"sim_velocity_{ID}.txt"))

    os.chdir(os.path.join(os.getcwd(), temp_dir))





    # Save results:
    os.chdir(os.path.join("..", ".."))
    sim_temp.to_csv(os.path.join(RESULTS_DIR,
                                 f"sim_temp_{ID}.txt"), sep=" ")

    # Nettoyage du répertoire scratch temp_{ID} (2026-07-22) : jamais supprimé
    # auparavant (ligne laissée en commentaire) - sur un grand balayage (512
    # combinaisons), l'accumulation de tous les temp_{ID} a rempli le disque
    # (100% plein, 0 octet libre) et fait échouer ~15% des simulations. Les
    # échecs sont restés silencieux car pool.starmap_async(...).join() (voir
    # plus bas) n'appelle jamais .get() sur le résultat, donc n'importe quelle
    # exception dans un worker est avalée sans remonter au process principal.
    # shutil.rmtree ici évite la récidive sur un balayage encore plus grand.
    shutil.rmtree(temp_dir, ignore_errors=True)
    return


def _n_worker_processes():
    """Nombre de processus à lancer en parallèle, plafonné à MAX_WORKERS."""
    try:
        n_available = len(os.sched_getaffinity(0))
    except AttributeError:
        n_available = os.cpu_count() or 1
    return max(1, min(MAX_WORKERS, n_available - 2))


if __name__ == "__main__":
    import os

    # RUN_GRID_SEARCH=False (config_lomos.py) : on saute les simulations et on
    # réutilise les sim_temp_*.txt déjà présents, pour juste relancer le misfit
    # et/ou les plots sur un calage déjà terminé.
    if RUN_GRID_SEARCH:
        os.makedirs(os.path.join(BASE_APP_DIR, "temp"), exist_ok=True)
        os.makedirs(os.path.join(BASE_APP_DIR, "results"), exist_ok=True)

        # Compile le binaire "ginette" s'il n'existe pas encore (2026-07-22) :
        # GINETTE_SENSI/ginette est dans .gitignore (binaire compilé, dépendant de
        # l'OS/l'architecture/la glibc - ne doit pas être versionné), donc absent
        # après un git clone frais sur une autre machine/cluster/compte. Sans ce
        # bloc, _copy_from_app("ginette") plus bas échoue avec un FileNotFoundError
        # sans indication de comment le résoudre. On compile une seule fois ici
        # (processus principal, avant le pool) plutôt que dans chaque worker, pour
        # éviter une recompilation redondante par simulation.
        _ginette_sensi_dir = os.path.join(BASE_APP_DIR, "GINETTE_SENSI")
        os.makedirs(_ginette_sensi_dir, exist_ok=True)
        if not os.path.isfile(os.path.join(_ginette_sensi_dir, "ginette")):
            print(f"Binaire 'ginette' absent de {_ginette_sensi_dir}, compilation...")
            _cwd_before_compile = os.getcwd()
            os.chdir(_ginette_sensi_dir)
            try:
                compile_ginette_src(REPO_ROOT)
            finally:
                os.chdir(_cwd_before_compile)

        if (Delete_sim=="True"):
            delete_sim_temp(RESULTS_DIR, os.path.join(BASE_APP_DIR, "temp"))
        # grid search in results/{POINT_NAME}/grid_search.csv
        grid = pd.read_csv(os.path.join(RESULTS_DIR, "grid_search.csv"), delimiter=";")

        # find which simulations are already done
        # if file sim_temp_ID.txt exists in results/{POINT_NAME}/, consider it done
        if (Delete_sim!="True"):
            done = sorted([int(f.split("_")[-1].split(".")[0])
                       for f in os.listdir(RESULTS_DIR)])[1:]
            remains = [i for i in grid.ID if i not in done]
        else:
            remains = grid.ID.tolist()

        # log_k/n/lam/c : lus depuis la grille s'ils sont calibrés (colonne présente
        # dans grid_search.csv), sinon repli sur REF_LOG_K/REF_N/REF_LAM/
        # REF_HEAT_CAPACITY (grid search à 1, 2, 3 ou 4 paramètres, voir
        # 1_define_grid_search.py). nb_zone et alt_thk sont lus PAR LIGNE depuis la
        # grille (colonnes écrites par 1_define_grid_search.py), pas depuis une
        # constante fixée ici.
        # log_k2/lam2/n2 (zone de surface) : lus depuis la grille s'ils existent
        # (grid search à 2 zones), sinon None (grid search homogène à 1 zone).
        params = [[r.ID,
                   getattr(r, 'log_k', REF_LOG_K),
                   getattr(r, 'n', REF_N),
                   getattr(r, 'lam', REF_LAM),
                   getattr(r, 'c', REF_HEAT_CAPACITY),
                   date_simul_bg, dt, nb_day, state, z_top, z_bottom, az, dz, dz_obs,
                   r.nb_zone, r.alt_thk,
                   getattr(r, 'log_k2', None), getattr(r, 'n2', None),
                   getattr(r, 'lam2', None), getattr(r, 'c', REF_HEAT_CAPACITY)]
                  for r in grid.itertuples() if r.ID in remains]
        to = time()
        with mp.Pool(processes=_n_worker_processes()) as pool:
            result = pool.starmap_async(run_ginette, params)
            pool.close()
            pool.join()
            # result.get() est INDISPENSABLE (pas juste pool.join()) : sans lui,
            # une exception levée dans un worker (ex: disque plein) est avalée
            # silencieusement et le script se termine comme si tout avait
            # réussi, avec des sim_temp_ID.txt manquants sans aucun message
            # d'erreur (découvert 2026-07-22 sur un balayage à 512 combinaisons).
            result.get()
        tf = time()
        print(f"Run time for {len(params)} simulations: {round(tf-to, 2)} s")
    else:
        print("RUN_GRID_SEARCH=False : simulations non relancées, réutilisation des résultats existants.")

    # Chaînage vers les étapes suivantes (flags dans config_lomos.py) : chaque
    # script ne déclenche que le suivant, pour ne pas dupliquer la logique
    # d'enchaînement (3_misfit.py se charge lui-même de RUN_PLOTS).
    if RUN_MISFIT:
        subprocess.run([sys.executable, str(SCRIPT_DIR / "3_misfit.py")], check=True)

