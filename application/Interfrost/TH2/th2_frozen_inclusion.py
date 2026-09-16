#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cas de test InterFrost TH2 "Frozen Inclusion" (Grenier et al. 2018, Adv. Water
Resour. 114, 196-218) : dégel d'une inclusion carrée gelée (-5°C, 0.333 m de
côté, centrée en (1, 0.5)) dans un domaine 3 m x 1 m à +5°C traversé par un
écoulement horizontal (gradient de charge 0, 3, 9 ou 15 %). Température +5°C
imposée en amont, flux conductif nul ailleurs.

Ginette modélise la moitié haute ou basse du domaine (symétrie par rapport à
l'axe y = 0.5) : 3 m x 0.5 m, plan de symétrie en haut.

Mesures de performance (PM) de l'intercomparaison, écrites par Ginette
(ytest=TH2) dans S_TH2.dat à chaque pas de sortie, pour le domaine complet :
  PM1 : température minimale du domaine (°C) et temps de seuil (Tmin = 0°C)
  PM2 : flux de chaleur net sortant du système (W, épaisseur 1 m)
  PM3 : volume d'eau liquide (m3)
Référence : résultats des 13 codes de l'intercomparaison (site InterFrost),
convertis dans reference/interfrost2018_th2_pm.csv.

Usage :
    python3 th2_frozen_inclusion.py              # les 4 gradients
    python3 th2_frozen_inclusion.py 3 9          # gradients 3 % et 9 %
    python3 th2_frozen_inclusion.py 3 --dt 5 --dx 0.00833 --tend 180000
Sorties dans results/ : th2_pm_GH<g>.csv, th2_field_GH<g>.csv (champ final),
th2_pm_GH<g>.png, th2_field_GH<g>.png, th2_summary.csv, ginette_GH<g>.log.
"""

# %% IMPORTS
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import matplotlib
if not os.environ.get("MPLBACKEND") and not os.environ.get("DISPLAY") and sys.platform.startswith("linux"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
SRC_PY = REPO_ROOT / "src" / "src_python"
if str(SRC_PY) not in sys.path:
    sys.path.insert(0, str(SRC_PY))

from Init_folders import compile_ginette_src
from Interfrost_2D import (read_pm_records, read_s_pts_2d, read_coordinates, load_reference,
                           envelope, threshold_time, compare_to_reference)

CASE = "TH2"
GINETTE_SENSI = HERE / "GINETTE_SENSI"
RESULTS_DIR = HERE / "results"
REFERENCE = HERE / "reference" / "interfrost2018_th2_pm.csv"
RESULTS_DIR.mkdir(exist_ok=True)

# %% GÉOMÉTRIE (fiche InterFrost TH2)
Lx, Ly = 3.0, 1.0            # domaine complet
Lcx, Lcy, Lsq = 1.0, 0.5, 0.333   # inclusion carrée
T_in, T_frozen = 5.0, -5.0

# %% DISCRÉTISATION (valeurs par défaut, modifiables en ligne de commande)
dx = float("%8.2E" % (1.0 / 120))   # m, maille carrée (360 x 60 mailles sur le demi-domaine) ; 8 caractères pour le format d8.0
dt = 30.0                    # s (10, 30 et 60 s donnent le meme temps de seuil a 0.6 % pres)
t_end = 180000.0             # s (50 h) : l'inclusion est entièrement dégelée bien avant, même à 0 %
itsortie = 100               # s, pas d'enregistrement des PM
solver = "ILU"               # BiCGSTAB + ILU(0) : les solveurs a preconditionneur diagonal stagnent sur le systeme de pression (contraste de permeabilite 1e-6)
GRADIENTS = [0, 3, 9, 15]    # %


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gradients", nargs="*", type=int, default=GRADIENTS, help="gradients de charge en %% (0, 3, 9, 15)")
    ap.add_argument("--dx", type=float, default=dx, help="taille de maille [m]")
    ap.add_argument("--dt", type=float, default=dt, help="pas de temps nominal [s]")
    ap.add_argument("--tend", type=float, default=t_end, help="durée simulée [s]")
    ap.add_argument("--itsortie", type=int, default=itsortie, help="pas d'enregistrement des PM [s]")
    ap.add_argument("--solver", default=solver, choices=["ILU", "BIC", "CGS", "BIS"], help="solveur lineaire (ysolv)")
    return ap.parse_args()


def fortran_dx(value):
    """dx est lu par Ginette en d8.0 (8 caractères, ex. 8.33E-03) : on utilise
    dans tout le script la valeur effectivement lue, pour que le maillage
    reconstruit ici (nc, nr, coordonnées) soit exactement celui de Ginette."""
    return float("%8.2E" % value)


def apply_args(a):
    global dx, dt, t_end, itsortie, solver
    suffix = ""
    if a.solver != solver: suffix += f"_{a.solver}"
    if fortran_dx(a.dx) != dx: suffix += f"_dx{a.dx*1e3:g}mm"
    if a.dt != dt: suffix += f"_dt{a.dt:g}s"
    if a.tend != t_end: suffix += f"_t{a.tend:g}s"
    dx, dt, t_end, itsortie, solver = fortran_dx(a.dx), a.dt, a.tend, a.itsortie, a.solver
    return suffix


def mesh():
    """Demi-domaine 3 m x 0.5 m, maille carrée dx. Ginette numérote les mailles
    ligne par ligne depuis le haut (plan de symétrie), x croissant."""
    nc = int(round(Lx / dx))
    nr = int(round(Ly / 2 / dx))
    return nc, nr, nc * nr


def write_inputs(grad_pct, steady=False):
    """Substitue les placeholders des templates E_*_bck.dat -> E_*.dat (largeurs
    de champ imposées par les formats de lecture Fortran).

    steady=True : écoulement en régime permanent (rp=0, thermique éteinte,
    perméabilité du champ gelé initial) pour produire le champ de pression
    initial, comme la fiche InterFrost le demande. steady=False : transitoire,
    qui repart de ce champ (ichi2=1) pour les cas avec écoulement."""
    nc, nr, nm = mesh()
    az = Ly / 2
    flow = grad_pct > 0
    h_left = grad_pct / 100 * Lx          # charge amont (m), aval = 0
    subs = {
        "E_parametre": {
            "[iec]": "1" if flow else "0",
            "[rp]": "1",
            # étape amont : thermique active mais pas de temps très courts, la
            # pression s'équilibre en ~1 s (β μ L²/k) alors que la thermique
            # bouge sur des heures -> champ de pression permanent, T ~ inchangée
            "[ith]": "1",
            "[dt]": ("%9.3E" % 0.1) if steady else ("%9.3E" % dt),   # d9.0
            "[nitt]": "%10d" % (10 if steady else int(round(t_end))),   # i10, unite = 1 s
            "[al]": "%8.2E" % Lx,               # d8.0
            "[az]": "%9.3E" % az,               # d9.0
            "[z_top]": "%9.3E" % az,            # d9.0
            "[nb_cell]": "%06d" % nm,           # i6
            "[nci]": "%05d" % nc,               # i5
            "[nri]": "%05d" % nr,               # i5
            "[dx]": "%8.2E" % dx,               # d8.0
            "[dz]": "%8.2E" % dx,               # d8.0
            "[itsortie]": "%08d" % itsortie,    # i8
            "[ysolv]": solver,                  # A3
            "[nm1]": "%05d" % 1, "[nm2]": "%05d" % 2,   # sans effet pour TH2
            "[nm5]": "%05d" % 5, "[nm6]": "%05d" % 6, "[nm7]": "%05d" % 7, "[nm8]": "%05d" % 8,
        },
        "E_cdt_aux_limites": {
            "[icl_lr]": "-2" if flow else "-1",   # charges imposées / flux nul
            "[h_left]": "%12.5E" % h_left,        # d12.0
        },
        "E_cdt_initiale": {
            # transitoire avec écoulement : repart du champ de pression permanent
            "[ichi2]": "1" if (flow and not steady) else "0",
        },
    }
    for bck, sub in subs.items():
        txt = (GINETTE_SENSI / f"{bck}_bck.dat").read_text()
        for k, v in sub.items():
            txt = txt.replace(k, v)
        leftover = [line.split()[0] for line in txt.splitlines() if "=[" in line.split()[0]]
        assert not leftover, f"placeholders non substitués dans {bck}_bck.dat : {leftover}"
        (GINETTE_SENSI / f"{bck}.dat").write_text(txt)
    (GINETTE_SENSI / "E_p_therm.dat").write_text((GINETTE_SENSI / "E_p_therm_bck.dat").read_text())
    assert nm < 1000000 and nc < 100000 and nr < 100000
    return nc, nr, nm


def extract_steady_pressure(nm):
    """Lit la pression (Pa, colonne 4) du dernier pas de S_pression_charge_
    temperature.dat et l'écrit dans E_pression_initiale.dat (une valeur par
    maille, dans l'ordre des mailles) pour le transitoire (ichi2=1)."""
    prof = pd.read_csv(GINETTE_SENSI / "S_pression_charge_temperature.dat",
                       sep=r"\s+", header=None).tail(nm)
    pr = prof.iloc[:, 3].to_numpy(dtype=float)
    np.savetxt(GINETTE_SENSI / "E_pression_initiale.dat", pr, fmt="%.6E")
    return pr


def run_ginette(tag):
    os.chdir(GINETTE_SENSI)
    compile_ginette_src(REPO_ROOT.as_posix(), flags=("-O2",))
    for pattern in ["S_*", "fort.*"]:
        for f in GINETTE_SENSI.glob(pattern):
            f.unlink()
    log = RESULTS_DIR / f"ginette_{tag}.log"
    t0 = time.time()
    with open(log, "w") as fl:
        rc = subprocess.call(["./ginette"], stdout=fl, stderr=subprocess.STDOUT)
    wall = time.time() - t0
    print(f"ginette terminé (code {rc}) en {wall/60:.1f} min, log : {log}")
    return rc, wall


def plot_pm(pm, grad, tag):
    """Trois mesures de performance : enveloppe des 13 codes de référence
    (gris), médiane (noir) et Ginette (rouge)."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    specs = [("PM1", "PM1_Tmin_C", "PM1 - température minimale [°C]"),
             ("PM2", "PM2_Jout_W", "PM2 - flux de chaleur net sortant [W]"),
             ("PM3", "PM3_Vwater_m3", "PM3 - volume d'eau liquide [m³]")]
    stats = {}
    for ax, (name, col, label) in zip(axes, specs):
        curves = load_reference(REFERENCE, name, grad / 100)
        for m, (tm, vm) in curves.items():
            ax.plot(tm / 3600, vm, color="0.75", lw=0.8, label="codes InterFrost 2018" if m == min(curves) else None)
        t = pm["t_s"].to_numpy()
        lo, med, hi, n = envelope(curves, t)
        ok = n >= 3
        ax.plot(t[ok] / 3600, med[ok], "k-", lw=1.2, label="médiane des codes")
        ax.plot(t / 3600, pm[col], color="tab:red", lw=1.6, label="Ginette")
        stats[name] = compare_to_reference(t, pm[col].to_numpy(), curves)
        ax.set_xlabel("Temps [h]")
        ax.set_ylabel(label)
        ax.set_xlim(0, t[-1] / 3600)
        ax.grid(alpha=0.3)
        ax.set_title(f"{name} - RMSE/médiane {stats[name]['rmse_median']:.3g}, "
                     f"{100*stats[name]['frac_in_envelope']:.0f} % du temps dans l'enveloppe", fontsize=9)
    axes[0].axhline(0, color="gray", ls=":")
    axes[0].legend(fontsize=8)
    fig.suptitle(f"InterFrost TH2 - gradient {grad} % - dx = {dx*1e3:.2f} mm, dt = {dt:g} s", fontsize=11)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"th2_pm_{tag}.png", dpi=150)
    plt.close(fig)
    return stats


def plot_field(tag, nc, nr, nm, t_field, temp, sw):
    """Champ de température final sur le domaine complet (le demi-domaine
    calculé est reflété par rapport au plan de symétrie)."""
    T = temp.reshape(nr, nc)          # ligne 1 = haut (plan de symétrie)
    full = np.vstack([T[::-1], T])    # bas -> haut : moitié miroir puis moitié calculée
    fig, ax = plt.subplots(figsize=(10, 3.8))
    im = ax.imshow(full, origin="lower", extent=[0, Lx, 0, Ly], cmap="coolwarm", vmin=-5, vmax=5, aspect="equal")
    ax.contour(np.linspace(dx / 2, Lx - dx / 2, nc), np.linspace(dx / 2, Ly - dx / 2, 2 * nr), full,
               levels=[0.0], colors="k", linewidths=0.8)
    rect = plt.Rectangle((Lcx - Lsq / 2, Lcy - Lsq / 2), Lsq, Lsq, fill=False, ec="k", ls="--", lw=0.8)
    ax.add_patch(rect)
    fig.colorbar(im, ax=ax, label="Température [°C]", shrink=0.9)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(f"TH2 {tag} - température à t = {t_field/3600:.1f} h (isotherme 0°C en noir, inclusion initiale en tirets)", fontsize=9)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"th2_field_{tag}.png", dpi=150)
    plt.close(fig)


def run_gradient(grad, suffix=""):
    tag = f"GH{grad}{suffix}"
    print(f"\n########## InterFrost TH2 - gradient de charge {grad} % ##########")
    nc, nr, nm = mesh()
    wall_steady = 0.0
    if grad > 0:
        # 1) écoulement permanent -> champ de pression initial (fiche InterFrost)
        write_inputs(grad, steady=True)
        print(f"maillage : {nc} x {nr} = {nm} mailles de {dx*1e3:.2f} mm ; "
              f"étape amont : équilibrage de l écoulement (pression permanente)")
        _, wall_steady = run_ginette(f"{tag}_permanent")
        extract_steady_pressure(nm)
    # 2) transitoire (repart du champ permanent si écoulement)
    nc, nr, nm = write_inputs(grad, steady=False)
    print(f"transitoire : dt = {dt:g} s, t_end = {t_end:.0f} s"
          + (" (pression initiale = régime permanent)" if grad > 0 else " (pas d'écoulement)"))
    rc, wall = run_ginette(tag)
    wall += wall_steady

    pm = read_pm_records(GINETTE_SENSI / "S_TH2.dat", CASE)
    pm.to_csv(RESULTS_DIR / f"th2_pm_{tag}.csv", index=False, float_format="%.6g")
    t_last = pm["t_s"].iloc[-1]
    t_thr = threshold_time(pm["t_s"], pm["PM1_Tmin_C"], 0.0, rising=True)
    ref_thr = [threshold_time(tm, vm, 0.0) for tm, vm in load_reference(REFERENCE, "PM1", grad / 100).values()]
    ref_thr = np.array([v for v in ref_thr if np.isfinite(v)])
    print(f"t final = {t_last:.0f} s ({t_last/3600:.1f} h) - cible {t_end:.0f} s")
    print(f"PM1 : Tmin(0) = {pm['PM1_Tmin_C'].iloc[0]:.3f} °C, Tmin(fin) = {pm['PM1_Tmin_C'].iloc[-1]:.3f} °C, "
          f"temps de seuil Tmin=0°C : Ginette {t_thr/3600:.2f} h, codes InterFrost "
          f"{ref_thr.min()/3600:.2f}-{ref_thr.max()/3600:.2f} h (médiane {np.median(ref_thr)/3600:.2f} h)")
    print(f"PM3 : V_eau(0) = {pm['PM3_Vwater_m3'].iloc[0]:.4f} m3, V_eau(fin) = {pm['PM3_Vwater_m3'].iloc[-1]:.4f} m3 "
          f"(dégel complet : {0.37*Lx*Ly:.4f} m3)")
    stats = plot_pm(pm, grad, tag)
    for name, st in stats.items():
        print(f"{name} : RMSE à la médiane des codes = {st['rmse_median']:.4g}, dans l'enveloppe "
              f"{100*st['frac_in_envelope']:.0f} % du temps, écart/largeur d'enveloppe = {st['rel_spread']:.2f}")

    t_field, pr, temp, sw = read_s_pts_2d(GINETTE_SENSI / "S_pts", nm)
    x, z = read_coordinates(GINETTE_SENSI / "E_coordonnee.dat", nm)
    pd.DataFrame({"x_m": x, "y_m": z, "T_C": temp, "Sw": sw, "p_Pa": pr}).to_csv(
        RESULTS_DIR / f"th2_field_{tag}.csv", index=False, float_format="%.6g")
    plot_field(tag, nc, nr, nm, t_field, temp, sw)

    row = dict(case=CASE, gradient_pct=grad, dx_m=dx, dt_s=dt, t_final_s=t_last, wall_min=wall / 60,
               rc=rc, t_threshold_h=t_thr / 3600, t_threshold_ref_median_h=np.median(ref_thr) / 3600,
               t_threshold_ref_min_h=ref_thr.min() / 3600, t_threshold_ref_max_h=ref_thr.max() / 3600,
               Vwater_final_m3=pm["PM3_Vwater_m3"].iloc[-1])
    for name, st in stats.items():
        row[f"{name}_rmse_median"] = st["rmse_median"]
        row[f"{name}_frac_in_envelope"] = st["frac_in_envelope"]
    return row


# %% MAIN
if __name__ == "__main__":
    args = parse_args()
    suffix = apply_args(args)
    rows = [run_gradient(g, suffix) for g in args.gradients]
    summary = pd.DataFrame(rows)
    summary.to_csv(RESULTS_DIR / f"th2_summary{suffix}.csv", index=False, float_format="%.6g")
    print("\n=== Récapitulatif ===")
    print(summary.to_string(index=False))
