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
from config_lomos import POINT_NAME, BOTTOM_SENSOR_TRIVIAL_INDEX, DOMAIN_LENGTH, MAX_WORKERS, RUN_PLOTS

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
    from Analytical_validation import (bulk_thermal_conductivity, volumetric_heat_capacity,
                                        darcy_velocity_from_head, CW_VOL)
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
    """Nombre de processus à lancer en parallèle, plafonné à MAX_WORKERS."""
    try:
        n_available = len(os.sched_getaffinity(0))
    except AttributeError:
        n_available = os.cpu_count() or 1
    return max(1, min(MAX_WORKERS, n_available - 2))

# %% MISFIT:
# in repertory results, for each simulation done, compute misfit with observed data
results = pd.read_csv(os.path.join(RESULTS_DIR,"grid_search.csv"), delimiter=";")
obs_data = pd.read_csv(os.path.join(RESULTS_DIR,"observed_data.txt"), delimiter=" ",
                       index_col=[0])



results["misfit_1"] = np.nan
results["misfit_2"] = np.nan
results["misfit_3"] = np.nan
results["mae_1"] = np.nan
results["mae_2"] = np.nan
results["mae_3"] = np.nan
results["pbias_1"] = np.nan
results["pbias_2"] = np.nan
results["pbias_3"] = np.nan
results["kge_1"] = np.nan
results["kge_2"] = np.nan
results["kge_3"] = np.nan
results["spin_up_days"] = np.nan

done = sorted([
    int(f.split("_")[-1].split(".")[0])
    for f in os.listdir(RESULTS_DIR)
    if f.startswith("sim_temp_") and f.endswith('.txt')
])

print("dans done il y a :", done)

# Misfit by simulations:
err = 1
results["err_misfit"] = err

# Index par ID pour un lookup rapide de la ligne grid_search dans
# spin_up_days_for_row (appelée une fois par simulation, dans chaque worker).
results_by_id = results.set_index("ID", drop=False)

# Densité fixe (E_p_therm.dat, rhosi=1180 - sédiment riche en matière
# organique, baissé depuis 2650/quartz pur le 2026-07-24 pour que c=2000
# J/kg/K donne une capacité volumique saturée cohérente avec la littérature
# streambed ; CASE('ZHZ') dans src/ginette_V2.f90 ne la spatialise jamais -
# voir 2_run_real_case.py). Repli si lam/n/c ne sont pas calibrés pour ce grid
# search (mêmes valeurs de référence que 2_run_real_case.py, meilleur point
# lomos230) : à garder synchronisé si ces défauts changent là-bas.
REF_DENSITY = 1180.0
REF_LAM = 1.57
REF_N = 0.51
REF_HEAT_CAPACITY = 2000.0


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
    """Jours de spin-up nécessaires pour la simulation i, dérivés du temps
    caractéristique de diffusion thermique sur la colonne (tau = L^2/Gamma,
    même raisonnement que l'analyse Stallman dans
    stallman_diffusivity.py : tau ~ 0.4^2/0.25e-6 ~ 640000s ~ 7-8j)
    au lieu d'une valeur fixe pour toutes les simulations - un matériau moins
    diffusif (lam bas, n haut) a besoin de plus de jours pour "oublier" un
    état initial imparfait, un matériau très diffusif de moins."""
    row = results_by_id.loc[i]
    lam = getattr(row, "lam", REF_LAM)
    n = getattr(row, "n", REF_N)
    c = getattr(row, "c", REF_HEAT_CAPACITY)
    k_bulk = bulk_thermal_conductivity(lam, n)
    Cv = volumetric_heat_capacity(n, REF_DENSITY, cs_specific=c)
    gamma = k_bulk / Cv
    tau_seconds = DOMAIN_LENGTH**2 / gamma
    return max(1, int(np.ceil(tau_seconds / 86400)))


def compute_misfit(i):
    """Calcule les misfits (L2, L1, PBIAS, KGE) pour la simulation i.

    Fonction de MODULE (pas une closure locale) : requis pour être picklable
    par multiprocessing. obs_data/err/RESULTS_DIR sont lus comme variables
    globales du module (fork sous Linux : héritées par copie des processus
    enfants, pas besoin de les repasser en argument).
    """
    sim = pd.read_csv(os.path.join(RESULTS_DIR, f"sim_temp_{i}.txt"),
                      delimiter=" ", index_col=[0])
    # Spin-up automatique (spin_up_days_for_row, dérivé de lam/n/c de CETTE
    # simulation) pendant lequel le modèle "oublie" un état initial imparfait :
    # on l'exclut de la comparaison plutôt que de biaiser le misfit avec la
    # transition. NB : remove_first_two_days_time_based renvoie des COPIES
    # filtrées, il faut donc bien récupérer sa valeur de retour (sinon le
    # filtrage n'a aucun effet).
    spin_up_days = spin_up_days_for_row(i)
    sim_f, obs_f = remove_first_two_days_time_based(sim, obs_data, days=spin_up_days)
    m1 = misfit_L2(obs_f.Temp1, sim_f.Temp1, err=err)
    m2 = misfit_L2(obs_f.Temp2, sim_f.Temp2, err=err)
    m3 = misfit_L2(obs_f.Temp3, sim_f.Temp3, err=err)
    a1 = misfit_L1(obs_f.Temp1, sim_f.Temp1, err=err)
    a2 = misfit_L1(obs_f.Temp2, sim_f.Temp2, err=err)
    a3 = misfit_L1(obs_f.Temp3, sim_f.Temp3, err=err)
    p1 = misfit_pbias(obs_f.Temp1, sim_f.Temp1)
    p2 = misfit_pbias(obs_f.Temp2, sim_f.Temp2)
    p3 = misfit_pbias(obs_f.Temp3, sim_f.Temp3)
    k1 = kge(obs_f.Temp1, sim_f.Temp1)
    k2 = kge(obs_f.Temp2, sim_f.Temp2)
    k3 = kge(obs_f.Temp3, sim_f.Temp3)
    return i, m1, m2, m3, a1, a2, a3, p1, p2, p3, k1, k2, k3, spin_up_days


# ==============================
# MAIN (calcul en parallèle)
# ==============================
if __name__ == "__main__":
    n_cores = _n_worker_processes()
    print(f"Utilisation de {n_cores} cœurs")

    with mp.Pool(processes=n_cores) as pool:
        results_list = list(tqdm(pool.imap(compute_misfit, done),
                                  total=len(done), desc="Compute misfit"))

    for i, m1, m2, m3, a1, a2, a3, p1, p2, p3, k1, k2, k3, spin_up_days in results_list:
        results.loc[results["ID"] == i, "misfit_1"] = m1
        results.loc[results["ID"] == i, "misfit_2"] = m2
        results.loc[results["ID"] == i, "misfit_3"] = m3
        results.loc[results["ID"] == i, "mae_1"] = a1
        results.loc[results["ID"] == i, "mae_2"] = a2
        results.loc[results["ID"] == i, "mae_3"] = a3
        results.loc[results["ID"] == i, "pbias_1"] = p1
        results.loc[results["ID"] == i, "pbias_2"] = p2
        results.loc[results["ID"] == i, "pbias_3"] = p3
        results.loc[results["ID"] == i, "kge_1"] = k1
        results.loc[results["ID"] == i, "kge_2"] = k2
        results.loc[results["ID"] == i, "kge_3"] = k3
        results.loc[results["ID"] == i, "spin_up_days"] = spin_up_days

    # Total misfit (quadratique) et likelihood
    results["misfit_tot"] = (results["misfit_1"] + results["misfit_2"] +
                             results["misfit_3"])

    results["Lik1"] = likelyhood(results["misfit_1"])
    results["Lik2"] = likelyhood(results["misfit_2"])
    results["Lik3"] = likelyhood(results["misfit_3"])

    # Total MAE : critère alternatif, moins sensible aux décalages de phase
    # (cf. Cognac & Ronayne 2023) - les deux sont disponibles pour comparaison,
    # ils ne désignent pas nécessairement la même combinaison comme "meilleure".
    results["mae_tot"] = results["mae_1"] + results["mae_2"] + results["mae_3"]

    # PBIAS et KGE moyennés (pas sommés comme misfit_tot/mae_tot) : ce sont des
    # métriques déjà normalisées (% et score borné à 1), pas des erreurs
    # cumulables entre profondeurs.
    results["pbias_mean"] = (results["pbias_1"] + results["pbias_2"] +
                             results["pbias_3"]) / 3
    results["kge_mean"] = (results["kge_1"] + results["kge_2"] +
                           results["kge_3"]) / 3

    # Si BOTTOM_SENSOR (config_lomos.py) tombe sur l'un des 3 capteurs de
    # comparaison (ex: config C), la maille correspondante compare à sa propre
    # condition limite - trivial, presque parfait par construction. Ça fausse
    # misfit_tot/kge_mean vers le haut si on ne l'exclut pas. misfit_fair est
    # le critère à regarder pour classer les simulations entre elles.
    if BOTTOM_SENSOR_TRIVIAL_INDEX is not None:
        _trivial_col = f"misfit_{BOTTOM_SENSOR_TRIVIAL_INDEX + 1}"
        results["misfit_fair"] = results["misfit_tot"] - results[_trivial_col]
    else:
        results["misfit_fair"] = results["misfit_tot"]

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
        best_id = results.loc[results["misfit_fair"].idxmin(), "ID"]
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

        thresh = results["misfit_fair"].min() * 1.05
        plateau = results[results["misfit_fair"] <= thresh]
        best = results.loc[results["misfit_fair"].idxmin()]

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
    if RUN_PLOTS:
        subprocess.run([sys.executable, os.path.join(BASE_APP_DIR, "4_plot_results.py")], check=True)