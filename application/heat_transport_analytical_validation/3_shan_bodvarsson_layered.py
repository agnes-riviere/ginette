#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cas 3/3 : validation de Ginette contre la solution PERMANENTE et
MULTI-COUCHES de Shan & Bodvarsson (2004) — extension de Bredehoeft &
Papadopulos (1965, cas 2/3, application/heat_transport_analytical_validation/2_...) à un
milieu à plusieurs couches de conductivité thermique distincte. Voir
src/src_python/Analytical_validation.py pour la théorie.

Principe : deux zones de conductivité thermique très différentes (via
NB_ZONE=2 de Ginette, comme dans le grid search 1D_Stream_aquifer_GridSearch),
même flux de Darcy constant sur tout le domaine (hypothèse du papier),
températures fixées en haut et en bas - comparaison du profil température
numérique (Ginette) au profil analytique à 2 couches.
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
from Analytical_validation import bulk_thermal_conductivity, shan_bodvarsson_layered_profile, CS_SPECIFIC

# %% PARAMÈTRES PHYSIQUES CONNUS - 2 ZONES
# Rappel (voir 1_define_grid_search.py) : dans generate_zone_parameters, la
# ZONE 2 est la couche DU HAUT (entre la surface et alt_thk), la ZONE 1 est
# la couche DU BAS (entre alt_thk et le fond). On garde le MÊME log_k dans
# les deux zones (même flux de Darcy constant sur tout le domaine, hypothèse
# du papier), seule la conductivité thermique change - contraste net (x4)
# pour rendre le changement de courbure à l'interface bien visible.
log_k = -13.0
poro_top, lam_top = 0.35, 4.0     # zone 2 (couche du haut)
poro_bot, lam_bot = 0.35, 1.0     # zone 1 (couche du bas)

T0 = 10.0   # température fixée en haut (z=0) [°C]
TB = 20.0   # température fixée en bas (z=L) [°C]

target_uz = 7.0e-7   # m/s, flux de Darcy cible, positif vers le bas (recharge)

# Géométrie du domaine : E_parametre_bck.dat a nmi=40/nri=40 CODÉS EN DUR -
# az/dz DOIT valoir exactement 40 (voir 2_bredehoeft_papadopulos_steady.py).
# Domaine réduit à 0.4 m (dz=0.01, résolution 2.5x plus fine que le 1er essai
# à L=1m/dz=0.025) pour tester si le résidu d'erreur concentré à l'interface
# vient de la résolution du maillage (voir discussion) - interface au milieu.
z_top = 0.0
z_bottom = -0.4
dz = 0.01  # 40 mailles (0.4/0.01)
alt_thk = -0.2
dz_obs = 0.1  # sans effet sur le régime permanent
az = abs(z_top - z_bottom)
L = az
b_top = abs(alt_thk - z_top)      # épaisseur de la couche du haut (zone 2)
b_bot = abs(z_bottom - alt_thk)   # épaisseur de la couche du bas (zone 1)

dt = 900
nb_day = 60  # domaine plus petit -> convergence plus rapide (temps caractéristique ~L^2)
date_simul_bg = pd.Timestamp("2000-01-01 00:00:00")

# %% PRÉDICTION ANALYTIQUE (Shan & Bodvarsson 2004, 2 couches)
k_bulk_top = bulk_thermal_conductivity(lam_top, poro_top)
k_bulk_bot = bulk_thermal_conductivity(lam_bot, poro_bot)
print("=== Prédiction analytique (Shan & Bodvarsson 2004, 2 couches) ===")
print(f"Couche du haut (zone 2) : épaisseur={b_top} m, k_bulk={k_bulk_top:.4f} W/m/K")
print(f"Couche du bas  (zone 1) : épaisseur={b_bot} m, k_bulk={k_bulk_bot:.4f} W/m/K")

# %% CHARGE HYDRAULIQUE NÉCESSAIRE POUR OBTENIR target_uz (même log_k partout)
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

# %% LANCEMENT DE GINETTE (régime permanent, state=0, NB_ZONE=2)
os.chdir(GINETTE_SENSI)
if not os.path.isfile("ginette"):
    print("Binaire 'ginette' absent, compilation...")
    compile_ginette_src(REPO_ROOT.as_posix())

z_obs = setup_ginette(dt, 0, nb_day, z_top, z_bottom, az, dz, date_simul_bg, dz_obs, amu=mu)
import shutil as _shutil
_shutil.copy("E_cdt_initiale_bck.dat", "E_cdt_initiale.dat")
initial_conditions(obs_temp, z_top, z_bottom, dz, z_obs)
boundary_conditions(obs_temp, dt)
# zone 1 = bas, zone 2 = haut (voir generate_zone_parameters/1_define_grid_search.py)
# 4e paramètre = capacité calorifique du solide (cpm), PAS la densité (voir
# 2_bredehoeft_papadopulos_steady.py pour le détail du changement Fortran
# 2026-07-22 sur CASE('ZHZ')).
generate_zone_parameters(z_bottom, dz, 2, alt_thk, log_k, poro_bot, lam_bot, CS_SPECIFIC,
                          REF_k2=log_k, REF_n2=poro_top, REF_l2=lam_top, REF_r2=CS_SPECIFIC)

import subprocess
subprocess.call(["./ginette"])

# %% FIGURE DU DOMAINE MODÉLISÉ (zones, paramètres) - lue DIRECTEMENT dans les
# fichiers écrits par Ginette (E_coordonnee.dat, E_zone.dat,
# E_zone_parameter.dat), PAS depuis les variables Python du script : montre
# ce que Ginette a réellement reçu, pas ce qu'on avait l'intention de lui
# donner (utile pour détecter une divergence entre les deux, par ex. un bug
# d'attribution de zone/maille comme celui trouvé et corrigé le 2026-07-22
# dans generate_zone_parameters).
coord = pd.read_csv("E_coordonnee.dat", sep=r"\s+", header=None, names=["id", "x", "z"])
zone_per_cell = pd.read_csv("E_zone.dat", header=None, names=["zone"])
coord["zone"] = zone_per_cell["zone"].to_numpy()
coord = coord.sort_values("z")

zone_params = pd.read_csv("E_zone_parameter.dat", sep=r"\s+", header=None,
                           names=["zone", "k", "n", "lam", "cpm"])

fig_dom, ax_dom = plt.subplots(figsize=(5.5, 7.5))
zone_colors = {1: "tab:blue", 2: "tab:orange"}
for zid in sorted(coord["zone"].unique()):
    sub = coord[coord["zone"] == zid]
    ax_dom.axhspan(sub["z"].min(), sub["z"].max(), color=zone_colors.get(zid, "gray"), alpha=0.25)
    p = zone_params.loc[zone_params["zone"] == zid].iloc[0]
    ax_dom.text(0.5, sub["z"].mean(),
                f"Zone {zid} ({sub['z'].max():.3f} à {sub['z'].min():.3f} m)\n"
                f"{len(sub)} mailles\n"
                f"k = {p['k']:.2e} m2\nn (poro) = {p['n']}\n"
                f"lam (solide) = {p['lam']} W/m/K\ncpm = {p['cpm']}",
                ha="center", va="center", fontsize=9,
                bbox=dict(boxstyle="round", fc="white", ec=zone_colors.get(zid, "gray")))

# Interface(s) entre zones : là où le zone change d'une maille à la suivante
z_sorted = coord["z"].to_numpy()
zone_sorted = coord["zone"].to_numpy()
for i in range(1, len(zone_sorted)):
    if zone_sorted[i] != zone_sorted[i - 1]:
        z_iface = (z_sorted[i] + z_sorted[i - 1]) / 2
        ax_dom.axhline(z_iface, color="black", linestyle="--", linewidth=1)
        ax_dom.text(0.02, z_iface, f"interface z≈{z_iface:.3f} m", fontsize=7,
                    ha="left", va="bottom", color="dimgray")

# Conditions limites : lues dans E_temp_t.dat / E_charge_t.dat (1ère ligne,
# constantes dans le temps pour ce cas permanent) plutôt que T0/TB/h_top du
# script.
temp_t = pd.read_csv("E_temp_t.dat", sep=r"\s+", header=None, names=["T_top", "T_bottom"])
charge_t = pd.read_csv("E_charge_t.dat", sep=r"\s+", header=None, names=["h_top", "h_bottom"])
T_top_file = temp_t["T_top"].iloc[0]
T_bottom_file = temp_t["T_bottom"].iloc[0]
h_top_file = charge_t["h_top"].iloc[0]
h_bottom_file = charge_t["h_bottom"].iloc[0]

ax_dom.annotate(f"T = {T_top_file}°C, h = {h_top_file:.4f} m", xy=(0.5, z_sorted.max()),
                 xytext=(0.5, z_sorted.max() + 0.03), ha="center", fontsize=9, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="tab:red"))
ax_dom.annotate(f"T = {T_bottom_file}°C, h = {h_bottom_file:.4f} m", xy=(0.5, z_sorted.min()),
                 xytext=(0.5, z_sorted.min() - 0.05), ha="center", fontsize=9, fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="tab:red"))

ax_dom.set_xlim(0, 1)
ax_dom.set_ylim(z_sorted.min() - 0.07, z_sorted.max() + 0.06)
ax_dom.set_xticks([])
ax_dom.set_ylabel("Altitude z [m]")
ax_dom.set_title("Domaine modélisé (lu dans les fichiers Ginette)\nShan & Bodvarsson (2004), 2 couches")
fig_dom.tight_layout()
fig_dom.savefig(RESULTS_DIR / "shan_bodvarsson_domaine.png", dpi=150)

# %% LECTURE DU PROFIL COMPLET (régime permanent -> S_pression_charge_temperature.dat,
# une ligne par maille - voir 2_bredehoeft_papadopulos_steady.py pour le détail du format
# et le piège des 10 colonnes, pas 6, du fichier)
profil = pd.read_csv("S_pression_charge_temperature.dat", sep=r"\s+", header=None,
                      names=["id", "x", "z", "pression", "charge", "T", "f1", "f2", "f3", "f4"])
last = profil.sort_values("z")

depth_ginette = -last["z"].to_numpy()  # conversion altitude -> profondeur positive vers le bas
T_ginette = last["T"].to_numpy()

# %% COMPARAISON NUMÉRIQUE / ANALYTIQUE
T_theorie = shan_bodvarsson_layered_profile(depth_ginette, target_uz, [b_top, b_bot],
                                             [k_bulk_top, k_bulk_bot], T0, TB)
err = T_ginette - T_theorie
rmse = np.sqrt(np.mean(err**2))
max_err = np.max(np.abs(err))
print(f"\n=== Comparaison Ginette vs théorie ===")
print(f"RMSE = {rmse:.4f} °C, erreur max = {max_err:.4f} °C (sur {T0}-{TB}°C, {TB-T0}°C d'écart total)")

comparison = pd.DataFrame({"depth_m": depth_ginette, "T_ginette": T_ginette, "T_theorie": T_theorie, "err": err})
comparison.to_csv(RESULTS_DIR / "shan_bodvarsson_comparison.csv", index=False)

# %% FIGURE
fig, ax = plt.subplots(figsize=(6, 8))
z_plot = np.linspace(0, L, 400)
ax.plot(shan_bodvarsson_layered_profile(z_plot, target_uz, [b_top, b_bot], [k_bulk_top, k_bulk_bot], T0, TB),
        -z_plot, label="Théorie (Shan & Bodvarsson)", color="black")
ax.scatter(T_ginette, -depth_ginette, color="tab:red", s=10, zorder=5, label="Ginette (régime permanent)")
ax.axhline(alt_thk, color="gray", linestyle=":", label=f"Interface des couches (z={alt_thk} m)")
ax.set_xlabel("Température [°C]")
ax.set_ylabel("Altitude [m]")
ax.legend()
ax.grid(alpha=0.3)
ax.set_title(f"Validation Ginette vs Shan & Bodvarsson (2004), 2 couches\nRMSE={rmse:.4f}°C")
fig.tight_layout()
fig.savefig(RESULTS_DIR / "shan_bodvarsson_comparison.png", dpi=150)
plt.show()

if max_err < 0.5:
    print("VALIDATION OK (erreur max < 0.5°C)")
else:
    print("ATTENTION : écart important - vérifier les paramètres/conventions.")
