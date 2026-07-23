#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cas 4 : reproduction de la Figure 2 de Kurylyk, Irvine, Carey, Briggs,
Werkema & Bonham (2017, Hydrological Processes 31, 2648-2661) —
comparaison de profils de température en régime PERMANENT pour des milieux
à 1, 2, 3 et 4 couches de conductivité thermique distincte, à 2 valeurs de
flux de Darcy (donc de nombre de Péclet thermique) par configuration.

Kurylyk et al. (2017) utilisent SUTRA (éléments finis) comme référence
numérique pour valider leur solution analytique multi-couches (extension de
Shan & Bodvarsson 2004, voir 3_shan_bodvarsson_layered.py et
src/src_python/Analytical_validation.py). On reproduit ici la même figure
en utilisant GINETTE (différences finies) à la place de SUTRA : test de
non-régression supplémentaire, plus exigeant que le cas 3 (jusqu'à 4
couches au lieu de 2), pour le fix matt() (moyenne harmonique de
conductivité thermique aux interfaces, 2026-07-23 — voir
project_matt_conductivity_interface_fix dans la mémoire du projet).

Paramètres de Table 1 de Kurylyk et al. (2017) : conductivités BULK et
épaisseurs RELATIVES des couches conservées EXACTEMENT. Le domaine total
est réduit de 100 m (papier) à 1 m — même logique que
2_bredehoeft_papadopulos_steady.py / 3_shan_bodvarsson_layered.py : sur
100 m, le temps caractéristique de diffusion (L²/gamma) se compterait en
siècles. Le flux de Darcy q est recalculé pour chaque scénario afin de
conserver EXACTEMENT les mêmes nombres de Péclet thermiques que le papier
(Pe = q*sum(bi/gamma_i), eq. 12 de Kurylyk et al. 2017, voir
Analytical_validation.thermal_peclet_number) : le profil T(z)/L a alors
une forme rigoureusement identique à celle de la Figure 2 originale, à une
échelle spatiale/temporelle différente.
"""

# %% IMPORTS
import sys
from pathlib import Path
import os
import shutil
import subprocess

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

from Direct_model import setup_ginette, initial_conditions, boundary_conditions
from Init_folders import compile_ginette_src
from Analytical_validation import shan_bodvarsson_layered_profile, thermal_peclet_number, CW_VOL
from Stallman_analysis import solid_thermal_conductivity

# %% SCÉNARIOS - Table 1 de Kurylyk et al. (2017)
# lam_bulk : conductivités BULK (W/m/K), b_frac : épaisseurs relatives
# (fraction du domaine total), ordonnées du HAUT vers le BAS - Pe : les 2
# nombres de Péclet thermiques cibles (papier : q=0.2 et 1 m/an sur 100 m).
SCENARIOS = {
    1: {"lam_bulk": [3.21], "b_frac": [1.0], "Pe": [0.83, 4.13]},
    2: {"lam_bulk": [1.61, 3.21], "b_frac": [0.5, 0.5], "Pe": [1.24, 6.20]},
    3: {"lam_bulk": [0.803, 3.21, 6.42], "b_frac": [0.5, 0.25, 0.25], "Pe": [1.96, 9.81]},
    4: {"lam_bulk": [0.803, 1.61, 3.21, 6.42], "b_frac": [0.25, 0.25, 0.25, 0.25], "Pe": [1.55, 7.75]},
}

poro = 0.35     # même porosité pour toutes les couches (seule lam varie, comme le papier)
log_k = -13.0   # perméabilité UNIFORME (même flux de Darcy constant sur tout le domaine, hypothèse du papier)

T0 = 10.0   # température fixée en haut (z=0) [°C] - cohérent avec l'axe de la Fig. 2 du papier (10-15°C)
TB = 15.0   # température fixée en bas (z=L) [°C]

L = 1.0       # m, domaine réduit (voir docstring)
dz = 0.01     # m -> 100 mailles (compatible avec des couches en 1/2/4 -> 100/50/25 mailles)
z_top, z_bottom = 0.0, -L
az = abs(z_top - z_bottom)
dz_obs = 0.1  # sans effet sur le régime permanent

dt = 900
nb_day = 90   # >> temps caractéristique de diffusion L²/gamma le plus lent des 4 scénarios (~60 j)
date_simul_bg = pd.Timestamp("2000-01-01 00:00:00")


def write_multilayer_zones(b, lam_bulk, poro_layers, log_k):
    """
    Écrit E_coordonnee.dat / E_zone.dat / E_zone_parameter.dat pour un
    domaine à n couches (n quelconque, PAS limité à 2 comme
    Direct_model.generate_zone_parameters) - même format et même
    convention d'ordre de maille (id=1 en surface) que cette dernière,
    généralisés. Zone i = couche i, numérotées du HAUT vers le BAS (ordre
    naturel de b/lam_bulk en entrée), CONTRAIREMENT à generate_zone_parameters
    (zone 2 = haut, zone 1 = bas) - sans importance ici puisque
    E_zone_parameter.dat est réécrit intégralement à chaque appel.
    """
    n = len(b)
    d = np.cumsum(b)  # profondeur cumulée du bas de chaque couche

    zvalues = np.sort(np.arange(z_bottom + dz / 2, dz / 2, dz))[::-1]  # centre des mailles, surface -> fond
    depth = -zvalues  # profondeur positive vers le bas
    zone = np.searchsorted(d, depth, side="right") + 1  # 1..n
    zone = np.clip(zone, 1, n)

    coord = pd.DataFrame({"id": np.arange(1, len(zvalues) + 1), "x": 0.5, "z": zvalues})
    coord.to_csv("E_coordonnee.dat", index=False, sep=" ", header=False)
    pd.Series(zone).to_csv("E_zone.dat", index=False, header=False)

    k = 10 ** log_k
    lines = []
    for i in range(n):
        lam_solid = solid_thermal_conductivity(lam_bulk[i], poro_layers[i])
        lines.append(f"{i + 1}\t{k:.2e}\t{poro_layers[i]:.2f} {lam_solid:.4f} 1000.00")
    with open("E_zone_parameter.dat", "w") as f:
        f.write("\n".join(lines) + "\n")


# %% BOUCLE SUR LES 4 SCÉNARIOS x 2 VALEURS DE PÉCLET
os.chdir(GINETTE_SENSI)
if not os.path.isfile("ginette"):
    print("Binaire 'ginette' absent, compilation...")
    compile_ginette_src(REPO_ROOT.as_posix())

results = {}  # (n_layers, pe_idx) -> dict(depth_ginette, T_ginette, T_theorie, q, Pe_actual, b, lam_bulk)

for n_layers, cfg in SCENARIOS.items():
    b = L * np.array(cfg["b_frac"])
    lam_bulk = np.array(cfg["lam_bulk"], dtype=float)
    poro_layers = [poro] * n_layers

    for pe_idx, Pe_target in enumerate(cfg["Pe"]):
        gamma = lam_bulk / CW_VOL
        q = Pe_target / np.sum(b / gamma)  # inverse de thermal_peclet_number (eq. 12 Kurylyk et al. 2017)

        T_ref_celsius = (T0 + TB) / 2
        mu = 2.414e-5 * 10 ** (247.8 / (T_ref_celsius + 133.15))
        K = 10 ** log_k * 1000.0 * 9.81 / mu
        dh = q * L / K
        h_top, h_bottom = dh, 0.0

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

        z_obs = setup_ginette(dt, 0, nb_day, z_top, z_bottom, az, dz, date_simul_bg, dz_obs, amu=mu)
        shutil.copy("E_cdt_initiale_bck.dat", "E_cdt_initiale.dat")
        initial_conditions(obs_temp, z_top, z_bottom, dz, z_obs)
        boundary_conditions(obs_temp, dt)
        write_multilayer_zones(b, lam_bulk, poro_layers, log_k)

        print(f"--- {n_layers} couche(s), Pe={Pe_target} (q={q * 86400 * 365:.3f} m/an) ---")
        subprocess.call(["./ginette"])

        profil = pd.read_csv("S_pression_charge_temperature.dat", sep=r"\s+", header=None,
                              names=["id", "x", "z", "pression", "charge", "T", "f1", "f2", "f3", "f4"])
        last = profil.sort_values("z")
        depth_ginette = -last["z"].to_numpy()
        T_ginette = last["T"].to_numpy()

        T_theorie = shan_bodvarsson_layered_profile(depth_ginette, q, b, lam_bulk, T0, TB)
        Pe_actual = thermal_peclet_number(q, b, lam_bulk)
        err = T_ginette - T_theorie
        rmse = np.sqrt(np.mean(err ** 2))
        print(f"    Pe visé={Pe_target:.2f}, Pe obtenu={Pe_actual:.2f} -> RMSE = {rmse:.4f} °C")

        results[(n_layers, pe_idx)] = dict(depth=depth_ginette, T_ginette=T_ginette, T_theorie=T_theorie,
                                            q=q, Pe=Pe_actual, b=b, lam_bulk=lam_bulk, rmse=rmse)

# %% FIGURE (mise en page identique à la Figure 2 de Kurylyk et al. 2017 :
# colonnes = 1/2/3/4 couches, lignes = Pe faible / Pe fort)
fig, axs = plt.subplots(2, 4, figsize=(14, 8), sharex=True, sharey=True)
for n_layers in SCENARIOS:
    for pe_idx in (0, 1):
        ax = axs[pe_idx, n_layers - 1]
        r = results[(n_layers, pe_idx)]
        z_plot = np.linspace(0, L, 400)
        T_plot = shan_bodvarsson_layered_profile(z_plot, r["q"], r["b"], r["lam_bulk"], T0, TB)
        ax.plot(T_plot, z_plot, color="black", lw=1.5, label="Analytique")
        ax.plot(r["T_ginette"], r["depth"], color="tab:red", marker="o", ms=2, lw=0.8, label="Ginette")
        ax.invert_yaxis()
        ax.grid(alpha=0.3)
        ax.set_title(f"{n_layers} couche(s), Pe={r['Pe']:.2f}\nRMSE={r['rmse']:.3f}°C", fontsize=10)
        if pe_idx == 1:
            ax.set_xlabel("Température [°C]")
        if n_layers == 1:
            ax.set_ylabel(f"Profondeur [m] (Pe {'faible' if pe_idx == 0 else 'fort'})")
axs[0, 0].legend(loc="lower right", fontsize=8)
fig.suptitle("Reproduction de la Figure 2 de Kurylyk et al. (2017) — Ginette vs solution analytique\n"
             "(domaine réduit 100 m -> 1 m, mêmes Pe que le papier)")
fig.tight_layout()
fig.savefig(RESULTS_DIR / "kurylyk_figure2_reproduction.png", dpi=150)
plt.show()

# %% RÉCAPITULATIF
print("\n=== Récapitulatif ===")
summary = pd.DataFrame([
    {"n_couches": n, "Pe_vise": SCENARIOS[n]["Pe"][pe], "Pe_obtenu": results[(n, pe)]["Pe"],
     "RMSE_C": results[(n, pe)]["rmse"]}
    for n in SCENARIOS for pe in (0, 1)
])
print(summary.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
summary.to_csv(RESULTS_DIR / "kurylyk_figure2_reproduction.csv", index=False)

max_rmse = summary["RMSE_C"].max()
if max_rmse < 0.5:
    print(f"\nVALIDATION OK (RMSE max = {max_rmse:.4f}°C < 0.5°C)")
else:
    print(f"\nATTENTION : RMSE max = {max_rmse:.4f}°C - vérifier les paramètres/nb_day.")