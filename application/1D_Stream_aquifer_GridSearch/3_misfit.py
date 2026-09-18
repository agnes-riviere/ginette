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
import sys
from time import time
import shutil
import subprocess
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

# Config partagée avec 0_boundary_conditions_real_case.py / 2_run_real_case.py
# via config_lomos.py (voir README.md) : point d'observation, config BC (A/C)
# et géométrie du domaine, pour ne plus avoir à les resynchroniser à la main.
from config_lomos import (POINT_NAME, CALIB_SENSORS, sensor_index, DOMAIN_LENGTH, MAX_WORKERS,
                          RUN_PLOTS, ERR_MISFIT, SPIN_UP_RUN_DAYS)

# Capteurs de calage (noms physiques, ex. ['Temp2', 'Temp3'] en config C) et
# indices des colonnes de results.txt (misfit_2, misfit_3, ...). Les deux
# capteurs de CL ne sont jamais comparés (voir CALIB_SENSORS, config_lomos.py).
SENSORS = list(CALIB_SENSORS)
IDX = [sensor_index(sn) for sn in SENSORS]

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

try:
    from Analytical_validation import darcy_velocity_from_head, CW_VOL
except Exception:
    print("Error importing Analytical_validation from src_python.")

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
def misfit_L2(obs, sim, err):
    """Erreur quadratique pondérée (chi², type RMSE). Sensible aux petits
    décalages de phase/amplitude entre signaux quasi-sinusoïdaux (cf. Cognac &
    Ronayne 2023, qui recommandent MAE pour cette raison sur un cas similaire)."""
    return np.sum(((sim - obs) / err)**2)


def misfit_L1(obs, sim, err):
    """Erreur absolue moyenne (MAE), pondérée par l'incertitude de mesure.
    Moins sensible qu'un critère quadratique aux petits décalages de phase
    entre le signal simulé et observé (cf. Cognac & Ronayne 2023, §2.3 :
    plusieurs fonctions objectif testées, MAE retenu pour cette robustesse)."""
    return np.sum(np.abs(sim - obs) / err)

def misfit_pbias(obs, sim):
    """Percent bias (PBIAS, Moriasi et al. 2007 / hydroGOF::pbias) : biais
    moyen entre simulé et observé, en %. Positif = surestimation, négatif =
    sous-estimation, 0 = pas de biais systématique. Ne capture PAS la forme
    du signal (une simulation constante à la moyenne observée peut donner un
    PBIAS proche de 0)."""
    return 100 * np.sum(sim - obs) / np.sum(obs)


def kge(obs, sim):
    """Kling-Gupta Efficiency (Gupta et al. 2009, J. Hydrology). KGE=1 :
    simulation parfaite ; décompose l'erreur en 3 composantes indépendantes :
    - r : corrélation linéaire (forme/timing du signal)
    - alpha = std(sim)/std(obs) : ratio de variabilité (amplitude)
    - beta = mean(sim)/mean(obs) : ratio de biais (niveau moyen)
    KGE = 1 - sqrt((r-1)^2 + (alpha-1)^2 + (beta-1)^2)."""
    r = np.corrcoef(obs, sim)[0, 1]
    alpha = np.std(sim) / np.std(obs)
    beta = np.mean(sim) / np.mean(obs)
    return 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)


def likelyhood(m):
    return np.exp(-0.5*m)


def log_likelyhood(m):
    return np.log(likelyhood(m))


def _n_worker_processes():
    """Nombre de processus à lancer en parallèle, plafonné à MAX_WORKERS
    (sauf si MAX_WORKERS="auto" : alors coeurs disponibles - 2, sans plafond)."""
    try:
        n_available = len(os.sched_getaffinity(0))
    except AttributeError:
        n_available = os.cpu_count() or 1
    if MAX_WORKERS == "auto":
        return max(1, n_available - 2)
    return max(1, min(MAX_WORKERS, n_available - 2))

# %% MISFIT:
# in repertory results, for each simulation done, compute misfit with observed data
results = pd.read_csv(os.path.join(RESULTS_DIR,"grid_search.csv"), delimiter=";")
obs_data = pd.read_csv(os.path.join(RESULTS_DIR,"observed_data.txt"), delimiter=" ",
                       index_col=[0])



for _metric in ("misfit", "mae", "pbias", "kge"):
    for _i in IDX:
        results[f"{_metric}_{_i}"] = np.nan
results["spin_up_days"] = np.nan

done = sorted([
    int(f.split("_")[-1].split(".")[0])
    for f in os.listdir(RESULTS_DIR)
    if f.startswith("sim_temp_") and f.endswith('.txt')
])

print("dans done il y a :", done)

# Misfit by simulations:
err = ERR_MISFIT  # voir config_lomos.py (ERR_MISFIT)
results["err_misfit"] = err

def mean_hydraulic_gradient():
    """Gradient de charge moyen (m/m) réellement appliqué comme CL, lu dans
    E_charge_t.dat plutôt que recalculé depuis le CSV brut - c'est la valeur
    que Ginette a effectivement vue (déjà correcte quelle que soit la config,
    y compris le x0.75 de la config C)."""
    path = os.path.join(BASE_APP_DIR, "GINETTE_SENSI", "E_charge_t.dat")
    charge = pd.read_csv(path, sep=" ", header=None, names=["h_top", "h_bottom"])
    return (charge["h_top"] - charge["h_bottom"]).mean() / DOMAIN_LENGTH


def ginette_velocity_for_row(i, spin_up_days):
    """Vitesse de Darcy (m/s) réellement calculée par Ginette pour la
    simulation i (résolution numérique complète de l'écoulement, pas juste
    k*gradient), lue dans sim_velocity_{i}.txt - copie de
    Sim_velocity_profil_t.dat faite par 2_run_real_case.py (colonnes: temps[s],
    z[m], vzm[m/s]). Moyennée sur la profondeur et sur le temps après spin-up.
    Retourne None si le fichier n'existe pas (grid search lancé avant l'ajout
    de cette sauvegarde, ou SAVE_VELOCITY_PROFILES=False) - 3_misfit.py se
    rabat alors sur le calcul Darcy analytique pour cette ligne.

    Signe : convention Ginette (z positif vers le haut, voir
    src/ginette_V2.f90 ~l.3712, v=-K/mu*(dP/dz+rho*g)) - vzm > 0 = flux vers
    le HAUT (exfiltration), vzm < 0 = vers le BAS (infiltration). C'est LA
    référence du projet ; q_analytic (darcy_velocity_from_head) est negé au
    point d'appel pour matcher cette même convention dans darcy_flux_m_s."""
    path = os.path.join(RESULTS_DIR, f"sim_velocity_{i}.txt")
    if not os.path.exists(path):
        return None
    vel = pd.read_csv(path, sep=r"\s+", header=None, names=["time", "z", "vzm"])
    vel_f = vel[vel["time"] >= spin_up_days * 86400]
    return vel_f["vzm"].mean() if not vel_f.empty else None


def spin_up_days_for_row(i):
    """Jours de spin-up exclus du misfit pour la simulation i.

    Le spin-up est ajouté AVANT DATE_SIMUL_BG (config_lomos.py : la simulation
    démarre SPIN_UP_RUN_DAYS jours plus tôt), donc on exclut exactement ces
    jours-là pour toutes les simulations : la fenêtre comparée est la même
    pour toute la grille (DATE_SIMUL_BG -> DATE_SIMUL_BG + NB_DAY). Le
    dimensionnement du spin-up (fixe ou "auto" = SPIN_UP_TAU_MULTIPLIER x tau
    du matériau le moins diffusif de la grille) vit dans config_lomos.py."""
    return SPIN_UP_RUN_DAYS


def compute_misfit(i):
    """Calcule les misfits (L2, L1, PBIAS, KGE) pour la simulation i.

    Fonction de MODULE (pas une closure locale) : requis pour être picklable
    par multiprocessing. obs_data/err/RESULTS_DIR sont lus comme variables
    globales du module (fork sous Linux : héritées par copie des processus
    enfants, pas besoin de les repasser en argument).
    """
    sim = pd.read_csv(os.path.join(RESULTS_DIR, f"sim_temp_{i}.txt"),
                      delimiter=" ", index_col=[0])
    # Spin-up (SPIN_UP_RUN_DAYS, config_lomos.py) pendant lequel le modèle
    # "oublie" un état initial imparfait : on l'exclut de la comparaison plutôt
    # que de biaiser le misfit avec la transition. NB :
    # remove_first_two_days_time_based renvoie des COPIES
    # filtrées, il faut donc bien récupérer sa valeur de retour (sinon le
    # filtrage n'a aucun effet).
    spin_up_days = spin_up_days_for_row(i)
    sim_f, obs_f = remove_first_two_days_time_based(sim, obs_data, days=spin_up_days)
    # Une valeur par capteur de calage ; sim_f et obs_f ont les mêmes noms de
    # colonnes (noms physiques des capteurs, ex. Temp2/Temp3 en config C).
    out = {}
    for sn, idx in zip(SENSORS, IDX):
        out[f"misfit_{idx}"] = misfit_L2(obs_f[sn], sim_f[sn], err=err)
        out[f"mae_{idx}"] = misfit_L1(obs_f[sn], sim_f[sn], err=err)
        out[f"pbias_{idx}"] = misfit_pbias(obs_f[sn], sim_f[sn])
        out[f"kge_{idx}"] = kge(obs_f[sn], sim_f[sn])
    out["spin_up_days"] = spin_up_days
    return i, out


# ==============================
# MAIN (calcul en parallèle)
# ==============================
if __name__ == "__main__":
    n_cores = _n_worker_processes()
    print(f"Utilisation de {n_cores} cœurs")

    with mp.Pool(processes=n_cores) as pool:
        results_list = list(tqdm(pool.imap(compute_misfit, done),
                                  total=len(done), desc="Compute misfit"))

    for i, out in results_list:
        for col, val in out.items():
            results.loc[results["ID"] == i, col] = val

    # Totaux sur les seuls capteurs de calage (les CL ne sont pas comparées, il
    # n'y a donc rien de trivial à retirer : misfit_tot est LE critère de
    # classement). Sommes pour misfit (L2) et mae (L1) : toutes les simulations
    # sont comparées sur la même fenêtre et le même nombre de points (spin-up
    # fixe avant DATE_SIMUL_BG), les sommes sont donc comparables entre elles.
    results["misfit_tot"] = sum(results[f"misfit_{_i}"] for _i in IDX)
    for _i in IDX:
        results[f"Lik{_i}"] = likelyhood(results[f"misfit_{_i}"])

    # Total MAE : critère alternatif, moins sensible aux décalages de phase
    # (cf. Cognac & Ronayne 2023) - les deux sont disponibles pour comparaison,
    # ils ne désignent pas nécessairement la même combinaison comme "meilleure".
    results["mae_tot"] = sum(results[f"mae_{_i}"] for _i in IDX)

    # PBIAS et KGE moyennés (pas sommés comme misfit_tot/mae_tot) : ce sont des
    # métriques déjà normalisées (% et score borné à 1), pas des erreurs
    # cumulables entre profondeurs.
    results["pbias_mean"] = sum(results[f"pbias_{_i}"] for _i in IDX) / len(IDX)
    results["kge_mean"] = sum(results[f"kge_{_i}"] for _i in IDX) / len(IDX)

    # Diagnostic d'identifiabilité (nombre de Péclet + largeur du plateau de
    # misfit sur log_k, voir README.md) : un bon misfit ne veut pas dire un
    # paramètre bien identifié, cf. Cucchi et al. 2021 (Pe<0.5 -> log_k pas
    # identifiable même si le point de misfit minimal a l'air net).
    try:
        grad = mean_hydraulic_gradient()
        # darcy_velocity_from_head(log_k, dh, L) divise dh par L en interne
        # (K*dh/L) - grad étant déjà dh/L (mean_hydraulic_gradient), il faut lui
        # repasser dh = grad*L, sinon L est divisé deux fois (q, donc Pe, était
        # surévalué d'un facteur 1/DOMAIN_LENGTH - trouvé en vérifiant le signe
        # de la vitesse Ginette, 2026-07-24).
        # darcy_velocity_from_head suit la convention "positif = infiltration"
        # (voir sa docstring/stallman_diffusivity.py) - opposée à la convention
        # Ginette (z positif vers le haut, positif = exfiltration, référence du
        # projet - voir ginette_velocity_for_row) : on negé ici pour que
        # darcy_flux_m_s reste cohérent quelle que soit sa source (analytique
        # ou Ginette).
        q_analytic = -darcy_velocity_from_head(results["log_k"].astype(float), grad * DOMAIN_LENGTH, DOMAIN_LENGTH)
        results["darcy_flux_m_s"] = q_analytic
        results["darcy_flux_source"] = "darcy_analytique"

        # Vitesse réellement calculée par Ginette (résolution numérique complète),
        # plus précise que la loi de Darcy analytique simplifiée, mais lue UNIQUEMENT
        # pour la meilleure simulation : sim_velocity_{ID}.txt fait ~9.5 Mo chacun,
        # les lire pour les 576 lignes (~5.4 Go, séquentiel) a fait planter la machine
        # (2026-07-24). Le reste des lignes garde le calcul analytique, déjà vectorisé.
        best_id = results.loc[results["misfit_tot"].idxmin(), "ID"]
        q_ginette_best = ginette_velocity_for_row(best_id, spin_up_days_for_row(best_id))
        if q_ginette_best is not None:
            results.loc[results["ID"] == best_id, "darcy_flux_m_s"] = q_ginette_best
            results.loc[results["ID"] == best_id, "darcy_flux_source"] = "ginette"

        # thermal_peclet_number() n'est pas vectorisable ligne à ligne (conçue pour
        # un système à N couches, un seul np.sum sans axis) - on applique sa formule
        # directement, équivalente pour une colonne homogène (1 couche) :
        # Pe = q*L/gamma = q*L*Cw_vol/lam.
        results["peclet"] = results["darcy_flux_m_s"] * DOMAIN_LENGTH * CW_VOL / results["lam"]

        print(f"\nVitesse Ginette lue pour la meilleure simulation uniquement (ID={int(best_id)}) "
              "- le reste utilise le calcul Darcy analytique.")

        thresh = results["misfit_tot"].min() * 1.05
        plateau = results[results["misfit_tot"] <= thresh]
        best = results.loc[results["misfit_tot"].idxmin()]

        print("\n=== Diagnostic d'identifiabilité ===")
        print(f"Gradient de charge moyen : {grad:.4f} m/m")
        print(f"Meilleur point : log_k={best['log_k']:.1f}, lam={best['lam']:.2f}, "
              f"n={best['n']:.2f}, Pe={best['peclet']:.4f}")
        print(f"Points dans le top 5% du misfit : {len(plateau)}/{len(results)}")
        print(f"  plage log_k : [{plateau['log_k'].min():.1f}, {plateau['log_k'].max():.1f}]")
        print(f"  plage lam   : [{plateau['lam'].min():.2f}, {plateau['lam'].max():.2f}]")
        print(f"  plage n     : [{plateau['n'].min():.2f}, {plateau['n'].max():.2f}]")
        # Le régime d'identifiabilité (Cucchi et al. 2021) dépend de la MAGNITUDE
        # du flux, pas de son sens - le signe de peclet indique la direction
        # (convention Ginette : + = exfiltration, - = infiltration), pas le régime.
        abs_pe = abs(best["peclet"])
        sens = "exfiltration" if best["peclet"] > 0 else "infiltration"
        print(f"Sens du flux (convention Ginette, +haut/-bas) : {sens}")
        if abs_pe < 0.5:
            print("  -> |Pe|<0.5 (régime conductif) : log_k borné seulement, pas identifiable.")
        elif abs_pe > 5:
            print("  -> |Pe|>5 (régime advectif) : log_k identifiable, lam/n non identifiables.")
        else:
            print("  -> 0.5<|Pe|<5 (transition) : log_k et diffusivité effective identifiables.")

        with open(os.path.join(RESULTS_DIR, "identifiability_diagnostic.txt"), "w") as f:
            f.write(f"gradient_m_m {grad:.6f}\n")
            f.write(f"best_log_k {best['log_k']:.2f}\n")
            f.write(f"best_lam {best['lam']:.2f}\n")
            f.write(f"best_n {best['n']:.2f}\n")
            f.write(f"best_peclet {best['peclet']:.4f}\n")
            f.write(f"plateau_n_points {len(plateau)}\n")
            f.write(f"plateau_log_k_min {plateau['log_k'].min():.2f}\n")
            f.write(f"plateau_log_k_max {plateau['log_k'].max():.2f}\n")
    except Exception as e:
        print(f"Diagnostic d'identifiabilité non calculé ({e}).")

    results.to_csv(os.path.join(RESULTS_DIR,"results.txt"), sep=" ")

    # Chaînage vers 4_plot_results.py (flag RUN_PLOTS dans config_lomos.py) -
    # ne se déclenche que si ce script est lancé seul ou par 2_run_real_case.py.
    # MPLBACKEND=Agg : 4_plot_results.py ne fait qu'enregistrer des PNG dans ce
    # contexte (pipeline batch) - le backend interactif natif macOS/Cocoa
    # ("macosx") plante en SIGSEGV quand il est sollicité depuis un subprocess
    # non-interactif comme celui-ci.
    if RUN_PLOTS:
        plot_env = dict(os.environ, MPLBACKEND="Agg")
        subprocess.run([sys.executable, os.path.join(BASE_APP_DIR, "4_plot_results.py")],
                       check=True, env=plot_env)