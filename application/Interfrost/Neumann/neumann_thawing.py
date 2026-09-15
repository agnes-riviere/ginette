#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validation de Ginette contre les solutions analytiques du cas test InterFrost
TH1 (dégel 1D d'un sol gelé, conduction +/- advection) - Kurylyk, McKenzie,
MacQuarrie & Voss (2014), "Analytical solutions for benchmarking cold regions
subsurface water flow and energy transport models: One-dimensional soil thaw
with conduction and advection", Advances in Water Resources 70, 172-184.

Trois benchmarks recommandés par le papier (Table S1 du supplément, copiée
dans reference/kurylyk2014_tableS1.csv) :

  1. Neumann      : pas d'écoulement, Ts=+5°C, Ti=-5°C (solution exacte à
                    deux phases, front franc) - "Run 15".
  2. Lunardini 98 : écoulement descendant v=10 m/an, Ts=+1°C, Ti=0°C
                    (solution quasi-stationnaire avec advection) - "Run 9".
  3. idem         : v=100 m/an - "Run 10".

Principe : colonne 1D initialement à Ti uniforme, température imposée en
surface Ts à t=0+ (palier), extrémité basse maintenue à Ti (loin du front,
approxime le milieu semi-infini). Ginette tourne avec ytest=TH1 (cas dédié :
force igel=2 et écrit la position des isothermes dans S_bound_permaf_1_t.dat)
et ymoycondtherm=NEUMA (conductivité bulk 2.619 / 1.839 W/m/K gelé / dégelé,
hardcodée dans ginette_V2.f90 = valeurs de la table InterFrost).

Les autres paramètres du benchmark sont reproduits UNIQUEMENT via les fichiers
d'entrée (aucune valeur hardcodée côté Fortran pour TH1, contrairement à THL) :
  C_dégelé = eps*rho_w*c_w + (1-eps)*rho_s*c_s  = 0.5*1000*4182 + 0.5*2500*889 = 3.201e6 J/m3/K
  C_gelé   = eps*rho_i*c_i + (1-eps)*rho_s*c_s  = 0.5*1000*2116 + 0.5*2500*889 = 2.169e6 J/m3/K
  L        = eps*rho_i*lat                       = 0.5*1000*334000            = 1.67e8  J/m3
Ginette calcule la chaleur latente en dégel avec rho_GLACE (rhoi*om*lat), alors
que le benchmark la définit avec rho_EAU (Swf*eps*rho_w*Lf) : on pose donc
rhoi=1000 et on ajuste cpice (2116 au lieu de 2300) pour garder C_gelé. La
solution analytique est appelée avec les valeurs effectivement vues par
Ginette (recalculées ci-dessous à partir des mêmes paramètres), pas avec les
valeurs arrondies de la table.

Usage :
    python3 neumann_thawing.py                 # benchmark 1 (Neumann) seul, 1 jour
    python3 neumann_thawing.py 1 2 3           # les trois benchmarks
    python3 neumann_thawing.py 1 --days 2 --dz 0.001 --dt 0.5 --tres -0.01   # sensibilité
Sorties dans results/ : neumann_front_<cas>.csv/.png, neumann_profile_<cas>.csv,
neumann_summary.csv, ginette_<cas>.log (avec un suffixe reprenant les options
si elles diffèrent des valeurs par défaut).
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
from Analytical_validation import (neumann_thaw_profile, kurylyk_advective_thaw_front,
                                   kurylyk_advective_thaw_profile)

GINETTE_SENSI = HERE / "GINETTE_SENSI"
RESULTS_DIR = HERE / "results"
REFERENCE = HERE / "reference" / "kurylyk2014_tableS1.csv"
RESULTS_DIR.mkdir(exist_ok=True)

# %% PARAMÈTRES PHYSIQUES (table InterFrost TH1 = Table 1/A1 de Kurylyk et al. 2014)
eps = 0.5          # porosité
rho_w, c_w = 1000.0, 4182.0   # eau
rho_s, c_s = 2500.0, 889.0    # solide
rho_i = 1000.0     # "glace" : = rho_w pour que rhoi*om*lat soit la chaleur latente du benchmark
c_i = 2116.0       # ajusté pour C_gelé = 2.169e6 J/m3/K (0.5*920*2300 avec la vraie densité de la glace)
Lf = 334000.0      # J/kg
k_u, k_f = 1.839, 2.619       # W/m/K, hardcodés dans ginette_V2.f90 (ymoycondtherm=NEUMA)
Tf = 0.0           # °C, température de fusion (tld/tlg)
YEAR = 365.25 * 86400.0

g_flow = 9.81      # m/s2, seulement pour les cas avec écoulement
k_perm = 1e-12     # m2, permeabilité intrinsèque (akx/akz du template) -> K = k*rho_w*g/mu
mu = 1e-3          # kg/m/s (amu du template)

C_u = eps * rho_w * c_w + (1 - eps) * rho_s * c_s
C_f = eps * rho_i * c_i + (1 - eps) * rho_s * c_s
L_vol = eps * rho_i * Lf

CASES = {
    1: dict(name="neumann", label="Benchmark 1 - Neumann (pas d'écoulement)",
            Ts=5.0, Ti=-5.0, v=0.0, ref_col="X_benchmark1_neumann_m"),
    2: dict(name="v10", label="Benchmark 2 - Lunardini 1998, v=10 m/an",
            Ts=1.0, Ti=0.0, v=10.0 / YEAR, ref_col="X_benchmark2_v10_m"),
    3: dict(name="v100", label="Benchmark 3 - Lunardini 1998, v=100 m/an",
            Ts=1.0, Ti=0.0, v=100.0 / YEAR, ref_col="X_benchmark3_v100_m"),
}

# %% DISCRÉTISATION (valeurs par défaut, modifiables en ligne de commande)
# Domaine : le benchmark est semi-infini. A t=nb_day, la perturbation thermique
# dans la zone gelée s'étend sur ~ 2*sqrt(alpha_f*t) = 0.64 m pour 1 jour ;
# la condition basse Dirichlet T=Ti doit être plus loin que ça.
# Fenêtre de phase T_res : le benchmark suppose un front franc (SUTRA :
# T_res=-0.0005°C avec dt=0.04-0.4 s). Ginette (capacité apparente "corde" +
# re-projection enthalpique, voir README) donne le même front pour
# T_res=-0.01 / -0.05 / -0.5°C et pour dt de 1 à 60 s ; seule dz compte
# (erreur ~dz/2 sur la position du front).
z_top = 0.0
z_bottom = -2.0
dz = 0.001      # m, comme la discrétisation SUTRA de Kurylyk et al. (2014)
dt = 10.0       # s
nb_day = 1
T_res = -0.01   # °C, bas de la fenêtre de changement de phase (tsg/tsd)
itsortie = 60   # s, pas d'enregistrement de S_bound_permaf_1_t.dat


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cases", nargs="*", type=int, default=[1], help="benchmarks à lancer (1, 2, 3)")
    ap.add_argument("--days", type=int, default=nb_day, help="durée simulée [j]")
    ap.add_argument("--dz", type=float, default=dz, help="taille de maille [m]")
    ap.add_argument("--dt", type=float, default=dt, help="pas de temps nominal [s]")
    ap.add_argument("--tres", type=float, default=T_res, help="bas de la fenêtre de phase tsg/tsd [°C]")
    ap.add_argument("--depth", type=float, default=abs(z_bottom), help="profondeur du domaine [m]")
    return ap.parse_args()


def apply_args(a):
    """Applique les options de la ligne de commande aux globales du module et
    renvoie un suffixe pour nommer les sorties (vide = valeurs par défaut)."""
    global dz, dt, nb_day, T_res, z_bottom, az, nb_cell, t_end
    suffix = ""
    if a.days != nb_day: suffix += f"_{a.days}d"
    if a.dz != dz: suffix += f"_dz{a.dz*1e3:g}mm"
    if a.dt != dt: suffix += f"_dt{a.dt:g}s"
    if a.tres != T_res: suffix += f"_tres{a.tres:g}"
    if a.depth != abs(z_bottom): suffix += f"_L{a.depth:g}m"
    dz, dt, nb_day, T_res, z_bottom = a.dz, a.dt, a.days, a.tres, -abs(a.depth)
    az = abs(z_top - z_bottom)
    nb_cell = int(round(az / dz))
    t_end = nb_day * 86400.0
    return suffix


az = abs(z_top - z_bottom)
nb_cell = int(round(az / dz))
t_end = nb_day * 86400.0


def write_inputs(case):
    """Substitue les placeholders des templates E_*_bck.dat -> E_*.dat.
    Les largeurs de champ sont imposées par les formats de lecture Fortran
    (lecture_parametre* dans ginette_V2.f90) : d8.0 -> 8 caractères, etc."""
    flow = case["v"] > 0
    # Ti "vu" par Ginette : le benchmark 2-3 part de Ti=Tf=0°C exactement, ce
    # qui pour Ginette (comme pour SUTRA, cf. note c de la table InterFrost)
    # serait un milieu déjà DÉGELÉ (sice=0 à T=tl). On part donc juste sous la
    # fenêtre de phase (Ti = T_res - 0.01), écart négligeable devant la chaleur
    # latente (C_f*0.06 / L ~ 1e-3).
    Ti_g = min(case["Ti"], T_res - 0.01)
    # Écoulement (benchmarks 2-3) : charges imposées en haut et en bas telles
    # que v = K*(h_top-h_bot)/az (loi de Darcy, K = k*rho_w*g/mu), écoulement
    # descendant. La charge initiale est le profil linéaire correspondant
    # (E_charge_initiale.dat), sinon le régime transitoire de pression (ss)
    # fausserait v pendant les premières ~ss*az^2/K secondes.
    K_hyd = k_perm * rho_w * g_flow / mu
    dh = case["v"] * az / K_hyd if flow else 0.0
    subs_param = {
        "[iec]": "1" if flow else "0",
        "[g]": "%8.2E" % (g_flow if flow else 0.0),   # d8.0
        "[dt]": "%9.3E" % dt,                 # d9.0
        "[nitt]": "%10d" % nb_day,            # i10 (unite = 86400 s)
        "[az]": "%9.3E" % az,                 # d9.0
        "[z_top]": "%9.2E" % z_top,           # d9.0
        "[z_bottom]": "%9.2E" % z_bottom,     # d9.0
        "[nb_cell]": "%06d" % nb_cell,        # i6 (nmi)
        "[nri]": "%05d" % nb_cell,            # i5 (nri)
        "[dz]": "%8.2E" % dz,                 # d8.0
        "[omp]": "%5.3f" % eps,               # f5.3
        "[itsortie]": "%08d" % itsortie,      # i8
    }
    subs_therm = {
        "[ithec]": "1" if flow else "0",
        "[rhoi]": "%8.2E" % rho_i,            # d8.0
        "[cpice]": "%10.4E" % c_i,            # d10.0
        "[lat]": "%10.4E" % Lf,               # d10.0
        "[tsd]": "%8.1E" % T_res,             # d8.0
        "[tld]": "%8.1E" % Tf,                # d8.0
    }
    subs_cl = {
        "[Ts]": "%12.5E" % case["Ts"],        # d12.0
        "[Ti]": "%12.5E" % Ti_g,              # d12.0
        "[icl_h]": "-2" if flow else "-1",    # i2 : charge imposée / flux nul
        "[icl_b]": "-2" if flow else "-1",
        "[h_top]": "%12.5E" % dh,             # d12.0, charge (m)
        "[h_bot]": "%12.5E" % 0.0,
    }
    subs_ci = {"[Ti]": "%8.1E" % Ti_g,        # d8.0
               "[ichi]": "1" if flow else "0"}
    if flow:
        z_cells = z_top - (np.arange(nb_cell) + 0.5) * dz
        h_cells = dh * (z_cells - z_bottom) / az
        np.savetxt(GINETTE_SENSI / "E_charge_initiale.dat", h_cells, fmt="%.8E")
    for bck, subs in [("E_parametre", subs_param), ("E_p_therm", subs_therm),
                      ("E_cdt_aux_limites", subs_cl), ("E_cdt_initiale", subs_ci)]:
        txt = (GINETTE_SENSI / f"{bck}_bck.dat").read_text()
        for k, v in subs.items():
            txt = txt.replace(k, v)
        leftover = [line.split()[0] for line in txt.splitlines() if "=[" in line.split()[0]]
        assert not leftover, f"placeholders non substitués dans {bck}_bck.dat : {leftover}"
        (GINETTE_SENSI / f"{bck}.dat").write_text(txt)
    # nri est lu en i5 : nb_cell doit tenir sur 5 chiffres
    assert nb_cell < 100000


def read_s_pts(path, nb_cell):
    """Lit S_pts (Fortran unformatted, sequential) : 1 enregistrement real(paso),
    puis nb_cell enregistrements (real(pr), real(temp), real(sw)). Chaque
    enregistrement est encadré par ses marqueurs de longueur (4 octets)."""
    raw = np.fromfile(path, dtype=np.uint8)
    assert len(raw) == 12 + nb_cell * 20, f"S_pts : {len(raw)} octets, attendu {12 + nb_cell*20}"
    paso = np.frombuffer(raw[4:8].tobytes(), dtype=np.float32)[0]
    body = np.frombuffer(raw[12:].tobytes(), dtype=np.float32).reshape(nb_cell, 5)
    return paso, body[:, 2].astype(float)


def analytical_front(case, t):
    if case["v"] == 0:
        return np.array([neumann_thaw_profile(0.0, ti, Ts=case["Ts"], Ti=case["Ti"], Tf=Tf,
                                              k_u=k_u, k_f=k_f, C_u=C_u, C_f=C_f, L=L_vol)[1]
                         for ti in np.atleast_1d(t)])
    return kurylyk_advective_thaw_front(t, case["v"], Ts=case["Ts"], Tf=Tf, k_u=k_u, C_u=C_u,
                                        L=L_vol, Cw_vol=rho_w * c_w)


def analytical_profile(case, depth, t):
    """Profil T(depth) analytique à l'instant t : Neumann (cas 1) ou
    quasi-stationnaire avec advection (cas 2-3)."""
    if case["v"] == 0:
        T, X, _ = neumann_thaw_profile(depth, t, Ts=case["Ts"], Ti=case["Ti"], Tf=Tf,
                                       k_u=k_u, k_f=k_f, C_u=C_u, C_f=C_f, L=L_vol)
    else:
        T, X = kurylyk_advective_thaw_profile(depth, t, case["v"], Ts=case["Ts"], Tf=Tf, k_u=k_u,
                                              C_u=C_u, L=L_vol, Cw_vol=rho_w * c_w)
    return T, X


def read_velocity(path, nb_cell):
    """Vitesse de Darcy verticale (m/s, >0 vers le haut dans Ginette) lue dans
    S_pression_charge_temperature.dat (colonnes vzm/vzp, E14.4)."""
    prof = pd.read_csv(path, sep=r"\s+", header=None).tail(nb_cell)
    return prof.iloc[:, 7].to_numpy(dtype=float), prof.iloc[:, 9].to_numpy(dtype=float)


def run_case(icase, reference, suffix=""):
    case = CASES[icase]
    tag = case["name"] + suffix
    print(f"\n########## {case['label']} ##########")
    os.chdir(GINETTE_SENSI)
    compile_ginette_src(REPO_ROOT.as_posix())
    write_inputs(case)
    for f in GINETTE_SENSI.glob("S_*"):
        f.unlink()
    for f in GINETTE_SENSI.glob("fort.*"):
        f.unlink()

    log = RESULTS_DIR / f"ginette_{tag}.log"
    t0 = time.time()
    with open(log, "w") as fl:
        rc = subprocess.call(["./ginette"], stdout=fl, stderr=subprocess.STDOUT)
    wall = time.time() - t0
    print(f"ginette terminé (code {rc}) en {wall/60:.1f} min, log : {log}")

    # %% FRONT DE DÉGEL X(t)
    # S_bound_permaf_1_t.dat (ytest=TH1) : t, zs(1,1), zs(1,2), zl(1,1), zl(1,2)
    # zs = altitude de l'isotherme tsd (bas de la fenêtre de phase), zl = celle
    # de tld (=Tf). -99 = isotherme absente. On prend l'isotherme Tf, que
    # Kurylyk et al. définissent comme la position du front ("shallowest node
    # below 0°C") ; avec T_res=-0.01°C les deux sont confondues à ~dz près.
    bound = pd.read_csv("S_bound_permaf_1_t.dat", sep=r"\s+", header=None,
                        names=["t", "zs1", "zs2", "zl1", "zl2"])
    bound = bound[bound["zl1"] > -90].copy()
    bound["X_ginette"] = z_top - bound["zl1"]
    bound["X_theorie"] = analytical_front(case, bound["t"].to_numpy())
    bound["err_m"] = bound["X_ginette"] - bound["X_theorie"]
    bound.to_csv(RESULTS_DIR / f"neumann_front_{tag}.csv", index=False)

    t_last = bound["t"].iloc[-1]
    X_last, X_th_last = bound["X_ginette"].iloc[-1], bound["X_theorie"].iloc[-1]
    late = bound[bound["t"] > 0.1 * t_last]   # les premiers instants sont dominés par dz
    rmse_X = np.sqrt(np.mean(late["err_m"]**2))
    rel_last = 100 * (X_last - X_th_last) / X_th_last
    print(f"t final atteint = {t_last:.0f} s ({t_last/86400:.3f} j) - cible {t_end:.0f} s")
    print(f"Front à t final : Ginette {X_last:.4f} m, théorie {X_th_last:.4f} m ({rel_last:+.2f} %)")
    print(f"RMSE sur X(t) pour t > 10 % de t final : {1e3*rmse_X:.2f} mm")

    # %% PROFIL DE TEMPÉRATURE AU DERNIER PAS DE SORTIE (S_pts, ytest=TH1)
    # S_pression_charge_temperature.dat est écrit en F14.2 (z=-0.00 sur un
    # maillage millimétrique, T à 0.01°C près) : on lit à la place S_pts, le
    # binaire réécrit à chaque pas de sortie avec T en simple précision. La
    # profondeur est reconstruite depuis l'indice de maille (maille 1 en haut).
    t_prof, T_g = read_s_pts(GINETTE_SENSI / "S_pts", nb_cell)
    depth = (np.arange(nb_cell) + 0.5) * dz
    assert abs(t_prof - t_last) < 1.0, f"S_pts (t={t_prof}) et S_bound_permaf_1_t.dat (t={t_last}) désynchronisés"
    summary = {"case": icase, "name": tag, "dz_m": dz, "dt_s": dt, "T_res_C": T_res, "t_final_s": t_last,
               "X_ginette_m": X_last, "X_theorie_m": X_th_last, "X_rel_err_pct": rel_last,
               "X_rmse_mm": 1e3 * rmse_X, "wall_min": wall / 60}
    T_th, _ = analytical_profile(case, depth, t_last)
    err = T_g - T_th
    rmse_T, max_T = np.sqrt(np.mean(err**2)), np.max(np.abs(err))
    print(f"Profil T à t={t_last:.0f} s : RMSE = {rmse_T:.4f} °C, erreur max = {max_T:.4f} °C")
    summary.update({"T_rmse_C": rmse_T, "T_max_err_C": max_T})
    prof = pd.DataFrame({"depth_m": depth, "T_ginette": T_g, "T_theorie": T_th, "err": err})
    if case["v"] > 0:
        vzm, vzp = read_velocity(GINETTE_SENSI / "S_pression_charge_temperature.dat", nb_cell)
        v_g = -np.median(vzm[1:-1])   # Ginette : vz > 0 vers le haut
        print(f"Vitesse de Darcy Ginette (médiane) = {v_g*YEAR:.3f} m/an, cible {case['v']*YEAR:.1f} m/an "
              f"(min {-vzm.max()*YEAR:.3f}, max {-vzm.min()*YEAR:.3f})")
        summary["v_ginette_m_per_yr"] = v_g * YEAR
        prof["vz_ginette"] = vzm
    prof.to_csv(RESULTS_DIR / f"neumann_profile_{tag}.csv", index=False)

    # %% FIGURE
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    t_plot = np.linspace(1.0, t_last, 400)
    ax1.plot(t_plot / 86400, analytical_front(case, t_plot), "k-", label="Théorie (Kurylyk et al. 2014)")
    ref = reference[reference["time_s"] <= t_last]
    ax1.plot(ref["time_days"], ref[case["ref_col"]], "o", mfc="none", color="gray", ms=5,
             label="Table S1, Kurylyk et al. 2014")
    ax1.plot(bound["t"] / 86400, bound["X_ginette"], "-", color="tab:red", lw=1.2, label="Ginette")
    ax1.set_xlabel("Temps [j]")
    ax1.set_ylabel("Profondeur du front de dégel X [m]")
    ax1.invert_yaxis()
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=8)
    ax1.set_title(f"{case['label']}\nX({t_last/86400:.2f} j) : Ginette {X_last:.4f} m / théorie {X_th_last:.4f} m ({rel_last:+.2f} %)",
                  fontsize=9)

    zmax = min(az, 4 * X_th_last)
    ax2.scatter(T_g, depth, s=6, color="tab:red", zorder=5, label="Ginette")
    d_plot = np.linspace(0, zmax, 500)
    T_plot, _ = analytical_profile(case, d_plot, t_last)
    ax2.plot(T_plot, d_plot, "k-", label="Théorie (Neumann)" if case["v"] == 0 else "Théorie (quasi-stationnaire, Kurylyk et al. 2014)")
    ax2.set_title(f"Profil à t={t_last/86400:.2f} j - RMSE={rmse_T:.3f} °C, max={max_T:.3f} °C", fontsize=9)
    ax2.axhline(X_th_last, color="gray", ls="--", label=f"Front théorique X={X_th_last:.4f} m")
    ax2.axvline(Tf, color="gray", ls=":")
    ax2.set_ylim(zmax, 0)
    ax2.set_xlabel("Température [°C]")
    ax2.set_ylabel("Profondeur [m]")
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"neumann_front_{tag}.png", dpi=150)
    plt.close(fig)
    return summary


# %% MAIN
if __name__ == "__main__":
    args = parse_args()
    suffix = apply_args(args)
    reference = pd.read_csv(REFERENCE, comment="#")
    rows = [run_case(i, reference, suffix) for i in args.cases]
    summary = pd.DataFrame(rows)
    summary.to_csv(RESULTS_DIR / f"neumann_summary{suffix}.csv", index=False)
    print("\n=== Récapitulatif ===")
    print(summary.to_string(index=False))
