#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validation de Ginette contre la solution analytique de Lunardini pour la
propagation d'un front de gel en régime TRANSITOIRE - cas "exact" (paramètres
non physiques, choisis pour coller aux hypothèses mathématiques de la
solution). Référence : V.J. Lunardini, "Freezing of soil with an unfrozen
water content and variable thermal properties", CRREL Report 82-2, 1988.
Voir src/src_python/Analytical_validation.py (lunardini_freezing_profile)
pour la théorie et le détail d'un signe du système d'équations.

Principe : colonne 1D initialement à T0=4°C uniforme, température imposée en
surface Ts=-6°C à t=0+ (palier), extrémité basse maintenue à T0 (loin du
front, approxime le milieu semi-infini de Lunardini). Pas d'écoulement
(charge uniforme) : test purement conductif avec changement de phase.

Ginette reproduit ce cas via 3 réglages dédiés de E_p_therm.dat :
- ymoycondtherm=LUNAR : conductivité bulk fixée à k1/k2/k3 (3.464352/2.941352/
  2.418352 W/m/K) selon T vs tsg (=Tm)/tlg (=Tf), zones gelée/mushy/dégelée.
- ytypsice=LINEA : teneur en glace linéaire entre tsg (sice=1-swressi) et
  tlg (sice=0) - modèle de Lunardini.
- ytest=THL : capacité volumique fixée à 690360 J/m3/K + terme de chaleur
  latente apparente (rhoi*om*lat*dSice/dT), reproduisant exactement le C
  constant à 3 zones et le alpha4 "effectif" de la solution analytique.
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
        cand_app = p / "application" / "Interfrost" / "Lunardi"
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
    BASE_APP_DIR = REPO_ROOT / "application" / "Interfrost" / "Lunardi"
if str(SRC_PY) not in sys.path:
    sys.path.insert(0, str(SRC_PY))

GINETTE_SENSI = BASE_APP_DIR / "GINETTE_SENSI"
RESULTS_DIR = BASE_APP_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

from Direct_model import setup_ginette, boundary_conditions, generate_zone_parameters
from Init_folders import compile_ginette_src
from Analytical_validation import lunardini_freezing_profile

# %% PARAMÈTRES DU CAS TEST (cas "exact" de Lunardini) - une colonne par valeur
# de Tm (solidus) testee ; tsg/tsd de E_p_therm_bck.dat portent un
# placeholder [Tm] substitue ci-dessous, donc pas besoin de synchroniser a
# la main comme avant.
Tm_list = [-1.0, -4.0]
T0 = 4.0     # °C, température initiale uniforme (= loin du front)
Ts = -6.0    # °C, température imposée en surface (Dirichlet, t>0)
poro = 0.336  # porosité (= theta du rapport, cf. gammad*xi0/rho_eau = 1680*0.2/1000)

# Géométrie : front à X2~0.25m (Tm=-1, 1 jour). Domaine porté à 3m (au lieu
# de 1m) - a 1m, la condition limite basse (Dirichlet, T=T0) contredisait
# la théorie elle-même (T3(x=1m) analytique ~= 2.93°C, pas 4°C : le milieu
# "semi-infini" n'est pas encore bien approxime a cette distance du front,
# erfc(beta*x/X2) pas encore negligeable). A x=3m, beta*x/X2 ~= 2.7,
# erfc(2.7) ~ 2e-4 : la theorie elle-meme est alors a T0 a moins de 0.01°C
# pres, la comparaison n'est plus polluee par cet artefact de bord.
z_top = 0.0
z_bottom = -3.0
dz = 0.005  # 600 mailles
az = abs(z_top - z_bottom)

dt = 5
nb_day = 1
date_simul_bg = pd.Timestamp("2000-01-01 00:00:00")
nb_cell = int(round(az / dz))

os.chdir(GINETTE_SENSI)
compile_ginette_src(REPO_ROOT.as_posix())

n_steps = int(nb_day * 86400 / dt) + 1
time_vector = np.arange(n_steps) * dt
t_actual = 86400.0

summary_rows = []
for Tm in Tm_list:
    print(f"\n########## Tm = {Tm}°C ##########")

    # %% LANCEMENT DE GINETTE (régime transitoire, state=1)
    # Palier de température (t=0+) et charge uniforme (pas d'écoulement) sur
    # toute la durée - PAS de gradient hydraulique, cas purement conductif.
    obs_temp = pd.DataFrame({
        "Time": time_vector,
        "T_top": np.full(n_steps, Ts),
        "T_bottom": np.full(n_steps, T0),
        "h_top": np.full(n_steps, 0.0),
        "h_bottom": np.full(n_steps, 0.0),
    })
    obs_temp.index = date_simul_bg + pd.to_timedelta(time_vector, unit="s")

    z_obs = setup_ginette(dt, 1, nb_day, z_top, z_bottom, az, dz, date_simul_bg, dz / 5, amu=1e-3)
    import shutil as _shutil
    _shutil.copy("E_cdt_initiale_bck.dat", "E_cdt_initiale.dat")

    # substitution du placeholder [Tm] (tsg/tsd) - setup_ginette() ne
    # substitue que [state] dans E_p_therm_bck.dat, le reste est a notre charge.
    with open("E_p_therm.dat") as f:
        therm = f.read()
    therm = therm.replace("[Tm]", "%8.3f" % Tm)
    with open("E_p_therm.dat", "w") as f:
        f.write(therm)

    # %% CONDITION INITIALE UNIFORME (T0 partout, charge nulle partout) - PAS
    # via initial_conditions() de Direct_model.py, qui interpole entre T_top et
    # T_bottom de la 1ere ligne de obs_temp (ici T_top=Ts, T_bottom=T0 : donnerait
    # un gradient au lieu d'un plateau uniforme a t=0). Meme format de sortie
    # (une valeur par ligne, maille du haut en premier - voir Direct_model.initial_conditions
    # et le fix du bug sep='\n' du 2026-09-10).
    with open("E_temperature_initiale.dat", "w") as f:
        pd.Series(np.full(nb_cell, T0)).to_csv(f, index=False, header=False)
    with open("E_charge_initiale.dat", "w") as f:
        pd.Series(np.full(nb_cell, 0.0)).to_csv(f, index=False, header=False)

    boundary_conditions(obs_temp, dt)
    # log_k tres faible : sans effet sur le resultat (charge uniforme -> vitesse
    # de Darcy nulle quel que soit K), garde juste generate_zone_parameters
    # valide. lam/cpm sans effet non plus (ymoycondtherm=LUNAR et ytest=THL
    # court-circuitent la conductivite/capacite generiques, voir E_p_therm_bck.dat).
    generate_zone_parameters(z_bottom, dz, 1, -0.5, -13.0, poro, 2.0, 840.0)

    import subprocess
    for f in GINETTE_SENSI.glob("S_pression_charge_temperature_day_*.dat"):
        f.unlink()
    subprocess.call(["./ginette"])

    # %% LECTURE DU PROFIL A ~1 JOUR (ytest=THL -> Sim_temperature_profil_t.dat
    # n'est JAMAIS ouvert : ce fichier n'existe que pour ytest in
    # ZHR/ZHZ/1DS/1DJ/ZND/ZNS/WAR (voir SELECT CASE(ytest) ~ligne 543 de
    # ginette_V2.f90). Le cas "THL" (dedie a Lunardini, voir E_p_therm_bck.dat)
    # ecrit a la place S_pression_charge_temperature_day_1.dat, un bloc de
    # nb_cell lignes reecrit a chaque pas de sortie - on prend le DERNIER bloc
    # (le plus proche de t=1 jour, cf. derniere ligne CONVERGE du run).
    sim = pd.read_csv("S_pression_charge_temperature_day_1.dat", sep=r"\s+", header=None,
                       names=["id", "x", "z", "pression", "charge", "Temp", "f1", "f2", "f3", "f4"])
    profil = sim.tail(nb_cell).sort_values("z")

    depth_ginette = -profil["z"].to_numpy()  # altitude -> profondeur positive vers le bas
    T_ginette = profil["Temp"].to_numpy()

    # %% COMPARAISON NUMÉRIQUE / ANALYTIQUE
    T_theorie, X1, X2 = lunardini_freezing_profile(depth_ginette, t_actual, Tm, T0=T0, Ts=Ts)
    err = T_ginette - T_theorie
    rmse = np.sqrt(np.mean(err**2))
    max_err = np.max(np.abs(err))
    print(f"=== Lunardini, Tm={Tm}°C, t={t_actual:.0f}s ===")
    print(f"Front gel/mushy X1={X1:.4f} m, front mushy/dégelé X2={X2:.4f} m")
    print(f"RMSE = {rmse:.4f} °C, erreur max = {max_err:.4f} °C (sur {Ts}-{T0}°C, {T0-Ts}°C d'écart total)")
    summary_rows.append({"Tm_C": Tm, "X1_m": X1, "X2_m": X2, "RMSE_C": rmse, "max_err_C": max_err})

    tag = "Tm_m%g" % abs(Tm) if Tm < 0 else "Tm_%g" % Tm
    tag = tag.replace(".", "_")
    comparison = pd.DataFrame({"depth_m": depth_ginette, "T_ginette": T_ginette, "T_theorie": T_theorie, "err": err})
    comparison.to_csv(RESULTS_DIR / f"lunardini_comparison_{tag}.csv", index=False)

    # %% FIGURE
    fig, ax = plt.subplots(figsize=(6, 8))
    z_plot = np.linspace(0.001, az, 400)
    T_plot, _, _ = lunardini_freezing_profile(z_plot, t_actual, Tm, T0=T0, Ts=Ts)
    ax.plot(T_plot, -z_plot, label="Théorie (Lunardini 1988)", color="black")
    ax.scatter(T_ginette, -depth_ginette, color="tab:red", s=8, zorder=5, label="Ginette")
    ax.axhline(-X1, color="gray", linestyle=":", label=f"Front gel/mushy (X1={X1:.3f} m)")
    ax.axhline(-X2, color="gray", linestyle="--", label=f"Front mushy/dégelé (X2={X2:.3f} m)")
    ax.set_xlabel("Température [°C]")
    ax.set_ylabel("Altitude [m]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title(f"Validation Ginette vs Lunardini (1988), Tm={Tm}°C\nRMSE={rmse:.4f}°C, t={t_actual/86400:.2f} j")
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"lunardini_comparison_{tag}.png", dpi=150)
    plt.close(fig)

    if max_err < 0.5:
        print("VALIDATION OK (erreur max < 0.5°C)")
    else:
        print("ATTENTION : écart important - vérifier les paramètres/conventions.")

# %% RÉCAPITULATIF
summary = pd.DataFrame(summary_rows)
summary.to_csv(RESULTS_DIR / "lunardini_summary.csv", index=False)
print("\n=== Récapitulatif ===")
print(summary.to_string(index=False))
