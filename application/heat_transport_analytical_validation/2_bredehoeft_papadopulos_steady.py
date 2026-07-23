#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cas 2/3 : validation de Ginette contre la solution PERMANENTE de
Bredehoeft & Papadopulos (1965) — profil température-profondeur en régime
permanent, milieu homogène, flux de Darcy vertical constant. Voir
src/src_python/Analytical_validation.py pour la théorie (cas particulier à
1 couche de Shan & Bodvarsson, voir 3_shan_bodvarsson_layered.py).

Principe : températures fixées en haut (T0) et en bas (TB) du domaine, flux
de Darcy constant imposé via une différence de charge - comparaison directe
du profil T(z) numérique (Ginette, régime permanent) au profil analytique.
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

from Direct_model import setup_ginette, initial_conditions, boundary_conditions, generate_zone_parameters
from Init_folders import compile_ginette_src
from Analytical_validation import bulk_thermal_conductivity, bredehoeft_papadopulos_profile, CW_VOL, CS_SPECIFIC

# %% PARAMÈTRES PHYSIQUES CONNUS
log_k = -13.0
lam = 2.0              # conductivité thermique de la fraction SOLIDE [W/m/K]
poro = 0.35

T0 = 10.0   # température fixée en haut (z=0) [°C]
TB = 20.0   # température fixée en bas (z=L) [°C]

target_uz = 7.0e-7   # m/s, flux de Darcy cible, positif vers le bas (recharge)

# Géométrie du domaine : E_parametre_bck.dat a nmi=40/nri=40 CODÉS EN DUR (pas
# de placeholder substitué par setup_ginette) - le nombre de mailles actives
# az/dz DOIT valoir exactement 40, sous peine de segfault (mailles Fortran
# dimensionnées pour 40, incohérence mémoire sinon). Domaine de 1 m choisi
# pour un temps caractéristique de diffusion (L^2/gamma) raisonnable - un
# domaine de plusieurs mètres (comme dans les exemples du papier) demanderait
# des centaines de jours de simulation pour converger en régime permanent.
z_top = 0.0
z_bottom = -1.0
dz = 0.025  # 40 mailles (1.0/0.025)
dz_obs = 0.2  # sans effet sur le régime permanent (pas de points fixes utilisés ici)
az = abs(z_top - z_bottom)
L = az

dt = 900
nb_day = 60  # ~2x le temps caractéristique de diffusion L^2/gamma pour ce domaine
date_simul_bg = pd.Timestamp("2000-01-01 00:00:00")

# %% PRÉDICTION ANALYTIQUE (Bredehoeft & Papadopulos 1965)
k_bulk = bulk_thermal_conductivity(lam, poro)
gamma = k_bulk / CW_VOL
Pe = target_uz * L / gamma
print("=== Prédiction analytique (Bredehoeft & Papadopulos 1965) ===")
print(f"k_bulk = {k_bulk:.4f} W/m/K, gamma = {gamma:.3e} m2/s, Pe = {Pe:.3f}")

# %% CHARGE HYDRAULIQUE NÉCESSAIRE POUR OBTENIR target_uz
T_ref_celsius = (T0 + TB) / 2
mu = 2.414e-5 * 10 ** (247.8 / (T_ref_celsius + 133.15))
K = 10**log_k * 1000.0 * 9.81 / mu
dh = target_uz * L / K
h_top, h_bottom = dh, 0.0
print(f"Charge imposée : h_top={h_top:.4f} m, h_bottom={h_bottom:.4f} m (K={K:.3e} m/s)")

# %% CONSTRUCTION DE LA CONDITION LIMITE (constantes dans le temps)
n_steps = int(nb_day * 86400 / dt) + 1
time_vector = np.arange(n_steps) * dt
obs_temp = pd.DataFrame({
    "Time": time_vector,
    "T_top": np.full(n_steps, T0),
    "T_bottom": np.full(n_steps, TB),
    "h_top": np.full(n_steps, h_top),
    "h_bottom": np.full(n_steps, h_bottom),
})
obs_temp.index = date_simul_bg + pd.to_timedelta(time_vector, unit="s")

# %% LANCEMENT DE GINETTE (régime permanent, state=0)
os.chdir(GINETTE_SENSI)
if not os.path.isfile("ginette"):
    print("Binaire 'ginette' absent, compilation...")
    compile_ginette_src(REPO_ROOT.as_posix())

z_obs = setup_ginette(dt, 0, nb_day, z_top, z_bottom, az, dz, date_simul_bg, dz_obs, amu=mu)
import shutil as _shutil
_shutil.copy("E_cdt_initiale_bck.dat", "E_cdt_initiale.dat")
initial_conditions(obs_temp, z_top, z_bottom, dz, z_obs)
boundary_conditions(obs_temp, dt)
# 4e paramètre = capacité calorifique du solide (cpm), PAS la densité :
# CASE('ZHZ') dans src/ginette_V2.f90 (2026-07-22) lit désormais cpmzone sur
# cette colonne de E_zone_parameter.dat ; rhos reste fixé à rhosi de
# E_p_therm.dat. CS_SPECIFIC (=1000) est sans effet sur ce cas (régime
# permanent : la prédiction analytique ne dépend pas de la capacité
# calorifique du solide), mais Ginette a quand même besoin d'une valeur.
generate_zone_parameters(z_bottom, dz, 1, -0.05, log_k, poro, lam, CS_SPECIFIC)

import subprocess
subprocess.call(["./ginette"])

# %% LECTURE DU PROFIL COMPLET (régime permanent : PAS de sortie temporelle
# Sim_temperature_maille*_t.dat / Sim_temperature_profil_t.dat - ces fichiers
# ne sont écrits qu'en régime TRANSITOIRE. En régime permanent (state=0),
# Ginette écrit directement le profil stabilisé, une ligne par maille, dans
# S_pression_charge_temperature.dat (colonnes : id, x, z, pression, charge,
# température, ... - voir Direct_model.reuse_end_in_initial qui en extrait
# les mêmes colonnes z/température pour amorcer un calcul transitoire).
# ATTENTION : ce fichier n'est PAS régénéré s'il existe déjà et que Ginette
# échoue silencieusement à le réécrire - vérifier sa date de modification
# après coup si un run précédent (transitoire ou avec un autre domaine) a pu
# laisser un fichier obsolète dans le même répertoire GINETTE_SENSI partagé.
profil = pd.read_csv("S_pression_charge_temperature.dat", sep=r"\s+", header=None,
                      names=["id", "x", "z", "pression", "charge", "T", "f1", "f2", "f3", "f4"])
last = profil.sort_values("z")

depth_ginette = -last["z"].to_numpy()  # conversion altitude -> profondeur positive vers le bas
T_ginette = last["T"].to_numpy()

# %% COMPARAISON NUMÉRIQUE / ANALYTIQUE
T_theorie = bredehoeft_papadopulos_profile(depth_ginette, target_uz, k_bulk, L, T0, TB)
err = T_ginette - T_theorie
rmse = np.sqrt(np.mean(err**2))
max_err = np.max(np.abs(err))
print(f"\n=== Comparaison Ginette vs théorie ===")
print(f"RMSE = {rmse:.4f} °C, erreur max = {max_err:.4f} °C (sur {T0}-{TB}°C, {TB-T0}°C d'écart total)")

comparison = pd.DataFrame({"depth_m": depth_ginette, "T_ginette": T_ginette, "T_theorie": T_theorie, "err": err})
comparison.to_csv(RESULTS_DIR / "bredehoeft_papadopulos_comparison.csv", index=False)

# %% FIGURE
fig, ax = plt.subplots(figsize=(6, 8))
z_plot = np.linspace(0, L, 200)
ax.plot(bredehoeft_papadopulos_profile(z_plot, target_uz, k_bulk, L, T0, TB), -z_plot,
        label="Théorie (Bredehoeft & Papadopulos)", color="black")
ax.scatter(T_ginette, -depth_ginette, color="tab:red", s=10, zorder=5, label="Ginette (régime permanent)")
ax.set_xlabel("Température [°C]")
ax.set_ylabel("Altitude [m]")
ax.legend()
ax.grid(alpha=0.3)
ax.set_title(f"Validation Ginette vs Bredehoeft & Papadopulos (1965)\nPe={Pe:.2f}, RMSE={rmse:.4f}°C")
fig.tight_layout()
fig.savefig(RESULTS_DIR / "bredehoeft_papadopulos_comparison.png", dpi=150)
plt.show()

if max_err < 0.5:
    print("VALIDATION OK (erreur max < 0.5°C)")
else:
    print("ATTENTION : écart important - vérifier les paramètres/conventions.")
