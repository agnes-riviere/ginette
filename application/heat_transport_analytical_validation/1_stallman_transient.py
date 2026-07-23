#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cas 1/3 : validation de Ginette contre la solution TRANSITOIRE de
Stallman (1965) — signal de température sinusoïdal en surface, milieu
homogène, écoulement de Darcy constant. Voir src/src_python/Analytical_validation.py
pour la théorie et OBS_point/notes_points_obs.txt (dossier voisin
1D_Stream_aquifer_GridSearch) pour le contexte du diagnostic qui a motivé
cette suite de validation.

Principe : on choisit des paramètres physiques CONNUS (log_k, lam, poro,
flux de Darcy cible), on en déduit l'amortissement/déphasage analytique
(a, b), on lance Ginette avec ces mêmes paramètres et une condition limite
de surface sinusoïdale pure (sans bruit), puis on ajuste un sinus sur la
sortie numérique à plusieurs profondeurs et on compare à la prédiction
théorique.
"""

# %% IMPORTS
import sys
from pathlib import Path
import os

import matplotlib
if not os.environ.get("MPLBACKEND") and not os.environ.get("DISPLAY") and sys.platform.startswith("linux"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Recherche les répertoires src/ et application/ du projet à partir de ce fichier
# (même mécanisme robuste que 1D_Stream_aquifer_GridSearch/0_boundary_conditions_real_case.py)
def find_project_paths(start_file=__file__):
    p = Path(start_file).resolve().parent
    repo_root = None
    src_py = None
    app_dir = None
    for _ in range(8):
        cand_src = p / "src" / "src_python"
        cand_app = p / "application" / "heat_transport_analytical_validation"
        if cand_src.exists():
            src_py = cand_src
        if cand_app.exists():
            app_dir = cand_app
        if src_py or app_dir:
            repo_root = p
            break
        p = p.parent
    return repo_root, src_py, app_dir


REPO_ROOT, SRC_PY, BASE_APP_DIR = find_project_paths()
if REPO_ROOT is None:
    REPO_ROOT = Path(__file__).resolve().parents[2]
if SRC_PY is None:
    SRC_PY = REPO_ROOT / "src" / "src_python"
if BASE_APP_DIR is None:
    BASE_APP_DIR = REPO_ROOT / "application" / "heat_transport_analytical_validation"
if str(SRC_PY) not in sys.path:
    sys.path.insert(0, str(SRC_PY))

GINETTE_SENSI = BASE_APP_DIR / "GINETTE_SENSI"
RESULTS_DIR = BASE_APP_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

from Direct_model import Direct_model_stallman
from Analytical_validation import (bulk_thermal_conductivity, volumetric_heat_capacity,
                                    stallman_forward, CW_VOL, CS_SPECIFIC)
from Stallman_analysis import fit_sinusoid

# %% PARAMÈTRES PHYSIQUES CONNUS (pilotent à la fois Ginette et la prédiction analytique)
log_k = -13.5          # log10(perméabilité intrinsèque k [m2]) - sable moyen à grossier
lam = 2.0              # conductivité thermique de la fraction SOLIDE [W/m/K]
poro = 0.35            # porosité [-]
REF_DENSITY = 2650.0   # densité du grain minéral [kg/m3] (quartz)

T_mean = 15.0          # température moyenne [°C]
A0 = 5.0               # amplitude du signal de surface [°C]
period_seconds = 86400.0  # cycle diurne

target_uz = 1.0e-6     # m/s, flux de Darcy cible, positif vers le bas (recharge)

# Géométrie du domaine : nmi/nri de E_parametre_bck.dat utilisent le
# placeholder [nb_cell] (plus de contrainte figée à 40 mailles). Domaine
# choisi assez profond pour que le signal diurne soit quasi totalement amorti
# au fond : 0.4 m - voir le test de sensibilité au domaine plus bas (0.4m
# donne un résultat équivalent à un domaine bien plus profond, avec un coût
# de calcul très inférieur).
z_top = 0.0
z_bottom = -0.4
dz = 0.01              # 40 mailles (0.4/0.01)
dz_obs = 0.1           # profondeurs d'observation : -0.1, -0.2, -0.3, -0.4 m
az = abs(z_top - z_bottom)

dt = 900               # s (15 min)
nb_day = 20            # 10 jours de mise en régime + 10 jours de fenêtre d'ajustement
SPINUP_DAYS = 10

date_simul_bg = pd.Timestamp("2000-01-01 00:00:00")

# %% PRÉDICTION ANALYTIQUE (théorie de Stallman - voir Analytical_validation.py)
k_bulk = bulk_thermal_conductivity(lam, poro)
Cv = volumetric_heat_capacity(poro, REF_DENSITY)
gamma = k_bulk / Cv
v_thermal = target_uz * CW_VOL / Cv
a_true, b_true = stallman_forward(v_thermal, gamma, period_seconds)

print("=== Prédiction analytique (Stallman 1965) ===")
print(f"k_bulk = {k_bulk:.4f} W/m/K, Cv = {Cv:.3e} J/m3/K, gamma = {gamma:.3e} m2/s")
print(f"v (vitesse frontale thermique) = {v_thermal:.3e} m/s")
print(f"a (amortissement) = {a_true:.4f} /m, b (déphasage) = {b_true:.4f} rad/m")

# %% VISCOSITÉ (référence physique du champ de pression interne SEULEMENT :
# depuis le passage au flux de Neumann imposé ci-dessous, amu n'affecte plus
# le flux réellement simulé - target_uz est imposé directement, indépendant
# de K/mu. Gardée cohérente avec T_mean par simple souci de réalisme physique
# du champ de pression, pas par nécessité pour la validation.)
mu = 2.414e-5 * 10 ** (247.8 / (T_mean + 133.15))

# %% CONDITION DE TEMPÉRATURE AU FOND DU DOMAINE
# T_bottom : imposer directement T_mean (signal totalement amorti) introduit
# un biais mesurable, le domaine (0.4 m, contraint à 40 mailles - voir plus
# haut) n'étant pas assez profond pour amortir totalement le cycle diurne
# (~5% de l'amplitude de surface subsiste à 0.4 m). On impose donc la valeur
# EXACTE prédite par la théorie à cette profondeur plutôt qu'une
# approximation "amorti à 100%" - cohérent puisque la solution de Stallman
# est justement ce qu'on cherche à valider par ailleurs. Testé (2026-07-23) :
# équivalent à un domaine profond (1m) avec T_bottom=T_mean simple (2.71%
# d'erreur contre 2.76% ici) - un domaine ENCORE plus profond (100m) est en
# revanche CONTRE-PRODUCTIF (32.8% d'erreur, le champ d'écoulement met trop
# longtemps à s'équilibrer sur un si grand domaine) - voir mémoire projet.
depth_bottom = abs(z_bottom - z_top)
amp_bottom = A0 * np.exp(-a_true * depth_bottom)
phase_bottom = b_true * depth_bottom

# %% LANCEMENT DE GINETTE (flux de Darcy imposé directement, condition de
# Neumann icl=-1 - voir Direct_model_stallman pour le détail et la
# convention de signe)
sim = Direct_model_stallman(GINETTE_SENSI, REPO_ROOT, dt, nb_day, z_top, z_bottom, dz, dz_obs,
                             date_simul_bg, log_k, poro, lam, CS_SPECIFIC, target_uz,
                             T_mean, A0, period_seconds, amp_bottom, phase_bottom, amu=mu)
depths = {1: 0.1, 2: 0.2, 3: 0.3}

# %% COMPARAISON NUMÉRIQUE / ANALYTIQUE (après exclusion du spin-up)
mask_fit = sim["Time"] >= SPINUP_DAYS * 86400
print(f"\n=== Comparaison Ginette vs théorie (après {SPINUP_DAYS} jours de spin-up) ===")
rows = []
for i, depth in depths.items():
    fit = fit_sinusoid(sim.loc[mask_fit, "Time"], sim.loc[mask_fit, f"Temp{i}"], period_seconds)
    amp_true = A0 * np.exp(-a_true * depth)
    phase_true = (b_true * depth) % (2 * np.pi)
    amp_err = 100 * (fit["amplitude"] - amp_true) / amp_true
    phase_err = fit["phase"] - phase_true
    rows.append({"depth_m": depth, "amp_ginette": fit["amplitude"], "amp_theorie": amp_true,
                 "err_amp_%": amp_err, "phase_ginette_rad": fit["phase"], "phase_theorie_rad": phase_true,
                 "err_phase_rad": phase_err})
comparison = pd.DataFrame(rows)
print(comparison.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
comparison.to_csv(RESULTS_DIR / "stallman_comparison.csv", index=False)

# %% FIGURE
fig, axs = plt.subplots(2, 1, figsize=(9, 8))
z_plot = np.linspace(0, 0.5, 200)
axs[0].plot(z_plot, A0 * np.exp(-a_true * z_plot), label="Amplitude théorique (Stallman)", color="black")
axs[0].scatter(comparison["depth_m"], comparison["amp_ginette"], color="tab:red", zorder=5, label="Ginette (ajusté)")
axs[0].set_xlabel("Profondeur [m]")
axs[0].set_ylabel("Amplitude diurne [°C]")
axs[0].legend()
axs[0].grid(alpha=0.3)
axs[0].set_title("Amortissement de l'amplitude avec la profondeur")

axs[1].plot(z_plot, (b_true * z_plot), label="Déphasage théorique (Stallman)", color="black")
axs[1].scatter(comparison["depth_m"], comparison["phase_ginette_rad"], color="tab:red", zorder=5, label="Ginette (ajusté)")
axs[1].set_xlabel("Profondeur [m]")
axs[1].set_ylabel("Déphasage [rad]")
axs[1].legend()
axs[1].grid(alpha=0.3)
axs[1].set_title("Déphasage avec la profondeur")

fig.suptitle("Validation Ginette vs Stallman (1965) — régime transitoire")
fig.tight_layout()
fig.savefig(RESULTS_DIR / "stallman_comparison.png", dpi=150)
plt.show()

max_amp_err = comparison["err_amp_%"].abs().max()
print(f"\nErreur relative maximale sur l'amplitude : {max_amp_err:.2f}%")
if max_amp_err < 10:
    print("VALIDATION OK (erreur < 10%)")
else:
    print("ATTENTION : écart important - vérifier les paramètres/conventions.")
