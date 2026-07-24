#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Estimation de la diffusivité thermique du sédiment et de la vitesse de Darcy
(infiltration/exfiltration) pour le point d'observation choisi (POINT_NAME dans
config_lomos.py, format LOMOS-mini), par la méthode de Stallman (1965) reprise
par Cheviron, Guérin, Tabbagh, Bendjoudi (2005) et décrite dans la note de
A. Tabbagh (26/06/2025) sur les profils du ru des Avenelles.

================================================================================
LE PRINCIPE, EN 4 ÉTAPES (à lire avant de dérouler le code)
================================================================================

Le ru chauffe/refroidit son lit chaque jour (cycle diurne). Cette onde de
chaleur se propage vers le bas : plus on descend, plus elle est AMORTIE
(amplitude plus faible) et plus elle est RETARDÉE (déphasage croissant). La
vitesse à laquelle amplitude et déphasage évoluent avec la profondeur dépend de
deux choses : (a) la diffusivité thermique du sédiment (à quelle vitesse la
chaleur diffuse), et (b) un éventuel écoulement d'eau vertical (qui "pousse" ou
"freine" la propagation de la chaleur selon son sens).

ÉTAPE 1 - Extraire l'onde diurne à chaque profondeur
  Pour chaque capteur (T_eau en surface, T1 à -10cm, ..., T4 à -40cm), on
  ajuste par moindres carrés un modèle T(t) = tendance_lente(t) + A*cos(wt) +
  B*sin(wt), avec w = 2*pi / 24h. On en tire, à chaque profondeur, une
  amplitude (A0) et une phase (à quel instant le maximum de chaleur passe-t-il).

ÉTAPE 2 - Regarder comment amplitude et phase évoluent avec la profondeur
  D'après la théorie, ln(amplitude) doit décroître LINÉAIREMENT avec la
  profondeur z, et la phase doit croître LINÉAIREMENT avec z. On trace donc
  deux droites (une régression chacune) : leurs pentes sont directement les
  grandeurs physiques a (amortissement) et b (déphasage) de la solution de
  Stallman T(z,t) = T0 + A0*exp(-a*z)*cos(wt - b*z).

ÉTAPE 3 - Un premier diagnostic : diffusivités apparentes
  Si on suppose qu'il n'y a AUCUN écoulement, a et b devraient être égaux et
  valoir tous les deux sqrt(w/(2*Gamma)). On calcule donc deux estimations
  indépendantes de la diffusivité Gamma, une à partir de a seul (Gamma_ampl),
  une à partir de b seul (Gamma_phase). Si elles sont proches : milieu homogène,
  sans écoulement notable, résultat fiable. Si elles divergent fortement : soit
  il y a un écoulement significatif, soit le milieu n'est pas homogène.

ÉTAPE 4 - Si Gamma_ampl != Gamma_phase : inverser Gamma et la vitesse ensemble
  On résout alors le système exact (a, b) -> (Gamma, v) sans supposer v=0, ce
  qui donne à la fois la diffusivité thermique ET la vitesse de Darcy (positive
  = infiltration vers le bas, négative = exfiltration vers le haut).

ÉTAPE 5 - De la vitesse de Darcy à la perméabilité (loi de Darcy)
  La colonne deltaP (H_corrige_24h) est la différence de charge hydraulique
  entre la rivière et un piézomètre installé à -40cm (positif = rivière plus
  haute que l'aquifère = infiltration attendue, cohérent avec le signe de uz).
  On a donc, sur la même colonne de sédiment (0 à -40cm) :

      uz = K * (dh/dz)     =>     K = uz / (dh/dz)    (loi de Darcy)

  K est la conductivité hydraulique [m/s]. Ce n'est PAS directement le
  paramètre calibré dans le grid search : 1_define_grid_search.py calibre
  log_k = log10(k), la PERMÉABILITÉ INTRINSÈQUE [m2], qui ne dépend que du
  milieu poreux (pas de la température de l'eau). On convertit donc :

      k = K * mu(T_eau) / (rho * g)

  où mu(T) est la viscosité de l'eau (dépend de la température du ru, mesurée
  par temperature_eau) et rho*g le poids volumique de l'eau.

ÉTAPE 5bis - Attention : Gamma_joint*Cv donne la conductivité du BULK, pas lam
  Gamma_joint * Cv donne la conductivité thermique du milieu SATURÉ (bulk, eau
  + solide mélangés) : ce n'est PAS lam. Le paramètre lam calibré dans le grid
  search (1_define_grid_search.py) est la conductivité thermique de la fraction
  SOLIDE seule (ks), une propriété du minéral qui ne dépend pas de la
  saturation en eau (confirmé : c'est bien ks que Ginette calibre en interne,
  voir alandas dans src/ginette_V2.f90). On convertit donc bulk -> solide, en
  utilisant la MÊME loi de mélange que Ginette (ymoycondtherm=ARITH dans
  E_p_therm.dat, moyenne arithmétique — pas le modèle "produit" de Cheviron et
  al. 2005/note de Tabbagh, que Ginette utiliserait seulement si configuré en
  ymoycondtherm=GEOME) :

      k_bulk = kw*n + ks*(1-n)     =>     ks = (k_bulk - kw*n) / (1 - n)

  où n est la porosité et kw = 0.598 W/m/K la conductivité thermique de l'eau.

ÉTAPE 6 - Sensibilité à la taille de fenêtre et à la porosité
  La porosité n est elle-même un paramètre calibré dans le grid search
  (n=[0.05, 0.5], pas fixé) : on balaie toute cette plage plutôt que d'en
  fixer une seule valeur. La taille de fenêtre (WINDOW_DAYS) est un choix
  méthodologique qui affecte fortement le résultat (voir plus bas) : on teste
  donc plusieurs tailles et on compare, plutôt que d'en choisir une au hasard.

================================================================================
Capteurs et profondeurs (dispositif LOMOS-mini, communes à lomos230/lomos231) :
  - temperature_eau (TempMolo) : 0 cm (référence de surface, eau du ru)
  - T1 : -10 cm, T2 : -20 cm, T3 : -30 cm, T4 : -40 cm (= profondeur du piézomètre deltaP)

Les points LOMOS-mini ne couvrent que quelques mois : on utilise donc le cycle
DIURNE (24h), pas le cycle annuel qui nécessiterait plusieurs années consécutives.

Voir src/src_python/Stallman_analysis.py pour les fonctions et la théorie
détaillée (formules, validation numérique de l'inversion).
================================================================================
"""

import os
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parents[2]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import numpy as np
import pandas as pd

# Backend non interactif si aucun affichage (appelé en pipeline depuis
# config_lomos.py, pas forcément depuis un terminal avec écran) - sinon
# plt.show() plante sans $DISPLAY au lieu de simplement ne rien afficher.
import matplotlib
if not os.environ.get("MPLBACKEND") and not os.environ.get("DISPLAY") and sys.platform.startswith("linux"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(project_root / "src" / "src_python"))
from Read_obs import read_lomos_csv
from Stallman_analysis import (fit_sinusoid, depth_profile_regression, darcy_law_K,
                                intrinsic_permeability, solid_thermal_conductivity)

APP_DIR = project_root / "application" / "1D_Stream_aquifer_GridSearch"
sys.path.insert(0, str(APP_DIR))
# POINT_NAME/SENSOR_DEPTHS/DELTAP_SENSOR partagés via config_lomos.py avec les
# autres scripts (0_/2_/3_/4_) - ce script analyse donc le même point et la
# même géométrie de capteurs que le reste du pipeline, pas toujours lomos230.
from config_lomos import POINT_NAME, SENSOR_DEPTHS, DELTAP_SENSOR

# =============================================================================
# PARAMÈTRES UTILISATEUR
# =============================================================================
OBS_DIR = APP_DIR / "OBS_point" / POINT_NAME
# Un sous-dossier par point (comme les autres scripts) : deux points analysés
# à la suite ne s'écrasent pas mutuellement.
RESULTS_DIR = APP_DIR / "results" / POINT_NAME
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Profondeurs des capteurs sous le fond du ru [m] (positif vers le bas)
DEPTHS = SENSOR_DEPTHS

# Distance entre la rivière (référence de charge) et le piézomètre où deltaP
# (H_corrige_24h) est mesuré.
DZ_HEAD = DEPTHS[DELTAP_SENSOR]

PERIOD = 86400.0          # période diurne [s] (24h)
TREND_ORDER = 3           # degré du polynôme de tendance lente (cubique, comme dans la note)

# Tailles de fenêtre testées [jours] : 30 vient de la note de Tabbagh, les
# autres explorent la sensibilité du résultat à ce choix (cf. ÉTAPE 6)
WINDOW_SIZES_TO_TEST = [7, 10, 15, 20, 25, 30]
# Nombre minimal de fenêtres exploitables pour qu'une taille soit considérée
# (en dessous, la statistique est trop pauvre pour être fiable)
MIN_WINDOWS_FOR_CANDIDATE = 5

# Porosité du sédiment saturé. C'est elle-même un paramètre calibré dans le grid
# search (n=[0.05, 0.5] dans 1_define_grid_search.py) : on ne la fixe donc pas
# à une seule valeur, mais on balaie toute sa plage pour voir comment les bornes
# suggérées pour lam et log_k se déplacent selon la porosité supposée. POROSITY
# sert seulement de valeur de référence pour le résumé détaillé fenêtre par fenêtre.
POROSITY = 0.5
POROSITY_RANGE = (0.05, 0.5)  # même plage que n dans le grid search

# Constantes alignées sur les valeurs par défaut internes de Ginette
# (GINETTE_SENSI/E_p_therm.dat), pour rester cohérent avec le modèle qui sera
# réellement calibré plutôt que d'utiliser des valeurs génériques de la littérature :
#   alandae=0598D-03  -> conductivité thermique eau = 0.598 W/m/K
#   cpe=4185D+00       -> chaleur spécifique eau = 4185 J/kg/K (x rho_eau=1000 -> volumique)
#   rhosi=1180D+00     -> masse volumique solide = 1180 kg/m3 (sédiment riche en matière
#                         organique, baissé depuis 2650/quartz pur le 2026-07-24)
#   cpm=2000           -> chaleur spécifique solide = 2000 J/kg/K (= c du grid search,
#                         config_lomos.py - ce champ E_p_therm.dat est mort côté Fortran
#                         pour ytest=ZHZ, écrasé par cpmzone, mais mis à jour pour
#                         cohérence avec ce qui est réellement simulé)
RHO_WATER = 1000.0
CW_WATER = 4185.0 * RHO_WATER   # chaleur volumique de l'eau [J/m3/K] (cpe * rho_eau)
RHO_SOLID = 1180.0
CS_SOLID = 2000.0 * RHO_SOLID   # chaleur volumique de la fraction solide [J/m3/K] (cpm * rhosi)
KW_WATER = 0.598                # conductivité thermique de l'eau [W/m/K] (alandae)

# Chaleur volumique du milieu saturé, dérivée de la porosité (mélange linéaire
# solide/eau, cf. Cheviron et al. 2005 éq. 7) : Cv = n*Cw + (1-n)*Cs
CV_SEDIMENT = POROSITY * CW_WATER + (1 - POROSITY) * CS_SOLID

# Seuil de désaccord relatif toléré entre Gamma_ampl et Gamma_phase pour qu'une
# fenêtre soit jugée "cohérente" (cf. ÉTAPE 3 ci-dessus)
COHERENCE_THRESHOLD = 0.5
# =============================================================================


print(__doc__)

# -----------------------------------------------------------------------------
# Lecture des données du point (voir Read_obs.read_lomos_csv : détecte et parse
# automatiquement ce format de fichier, quel que soit le format de date utilisé).
# On essaie chaque CSV du dossier du point jusqu'à en trouver un au format LOMOS-mini
# (même logique que process_obs_data dans 0_boundary_conditions_real_case.py) :
# le nom du fichier varie d'un point à l'autre (lomos230.csv, donnees_plot_lomos231.csv, ...).
# -----------------------------------------------------------------------------
df = None
csv_used = None
for csv_path in sorted(OBS_DIR.glob("*.csv")):
    df = read_lomos_csv(csv_path)
    if df is not None:
        csv_used = csv_path
        break
if df is None:
    raise ValueError(f"Aucun CSV au format LOMOS-mini trouvé dans {OBS_DIR} pour le point {POINT_NAME}.")
print(f"Fichier utilisé : {csv_used}")

# deltaP (H_corrige_24h) est enregistré en CENTIMÈTRES, pas en mètres (cf.
# 0_boundary_conditions_real_case.py, coef=0.01) : converti ici en mètres pour
# que le gradient dh/dz (ÉTAPE 5) et la loi de Darcy restent physiquement cohérents.
df['deltaP'] = df['deltaP'] * 0.01

df = df.set_index('dates').sort_index()
t0 = df.index.min()
time_s = (df.index - t0).total_seconds().to_numpy()
sensor_cols = list(DEPTHS.keys())

print(f"Période disponible : {df.index.min()} -> {df.index.max()} "
      f"({(df.index.max() - df.index.min()).days} jours)\n")


def analyze_window_size(window_days, verbose=True):
    """
    Découpe la série en fenêtres non chevauchantes de `window_days` jours et
    applique les ÉTAPES 1 à 5 sur chacune. Retourne (summary, window_data) :
    - summary : un DataFrame avec une ligne par fenêtre (résultats numériques).
    - window_data : dict {index_fenêtre: {'label','mask','fits','reg'}}, utile
      pour retracer les graphiques pédagogiques a posteriori.
    """
    window_seconds = window_days * 86400
    n_windows = int(time_s.max() // window_seconds)

    rows = []
    wdata = {}

    for w_idx in range(n_windows):
        t_start = w_idx * window_seconds
        t_end = t_start + window_seconds
        mask = (time_s >= t_start) & (time_s < t_end)
        if mask.sum() < 2 * window_days:  # sécurité : pas assez de points dans la fenêtre
            continue

        window_label = f"{(t0 + pd.Timedelta(seconds=t_start)).date()} -> {(t0 + pd.Timedelta(seconds=t_end)).date()}"

        # ÉTAPE 1 : un ajustement sinusoïdal indépendant par profondeur
        fits = {col: fit_sinusoid(time_s[mask], df[col].to_numpy()[mask], PERIOD, trend_order=TREND_ORDER)
                for col in sensor_cols}

        if verbose:
            print(f"=== Fenêtre {w_idx} : {window_label} ===")
            print("Étape 1 - amplitude et phase ajustées à chaque profondeur :")
            for name in sorted(DEPTHS, key=lambda n: DEPTHS[n]):
                retard_h = (fits[name]['phase'] % (2 * np.pi)) / (2 * np.pi) * 24
                print(f"  {name:9s} (z={DEPTHS[name]*100:4.0f} cm) : "
                      f"amplitude={fits[name]['amplitude']:.4f} °C, retard de phase={retard_h:.2f} h")

        # ÉTAPES 2 à 4 : régressions ln(amplitude) vs z et phase vs z, diffusivités
        # apparentes, puis inversion jointe (Gamma, vitesse de Darcy)
        reg = depth_profile_regression(fits, DEPTHS, PERIOD, Cv=CV_SEDIMENT, Cw=CW_WATER)

        if verbose:
            print("Étape 2 - régressions sur les 5 profondeurs :"
                  f" pente amplitude a={reg['a']:.3f} 1/m (R²={reg['r2_amp']:.3f}),"
                  f" pente phase b={reg['b']:.3f} rad/m (R²={reg['r2_phase']:.3f})")
            print(f"Étape 3 - diffusivités apparentes : Gamma_ampl={reg['gamma_ampl_1e-6']:.3f}e-6 m²/s, "
                  f"Gamma_phase={reg['gamma_phase_1e-6']:.3f}e-6 m²/s "
                  f"(désaccord relatif={reg['coherence_rel_diff']:.1%})")
            print(f"Étape 4 - inversion jointe : Gamma={reg['gamma_joint_1e-6']:.3f}e-6 m²/s, "
                  f"vitesse de Darcy uz={reg['uz_darcy_mm_yr']:.0f} mm/an "
                  f"({'infiltration' if reg['uz_darcy_mm_yr'] > 0 else 'exfiltration'})")

        # ÉTAPE 5 : gradient de charge (deltaP) sur la même fenêtre -> conductivité
        # hydraulique K, puis perméabilité intrinsèque k via la viscosité de l'eau
        dh_window = df['deltaP'].to_numpy()[mask]
        twater_window = df['TempMolo'].to_numpy()[mask]
        mean_dh = float(np.nanmean(dh_window))
        mean_twater = float(np.nanmean(twater_window))
        gradient_stable = bool(np.all(np.sign(dh_window) == np.sign(mean_dh))) if mean_dh != 0 else False
        # dh>0 (rivière plus haute) doit correspondre à une infiltration (uz>0), et
        # inversement : si les deux méthodes indépendantes (thermique et charge) ne
        # sont pas d'accord sur le SENS de l'écoulement, K n'a pas de sens physique
        # (un simple abs() masquerait ce désaccord au lieu de le signaler).
        signs_agree = bool(np.isfinite(reg['uz_darcy_m_s']) and np.sign(mean_dh) == np.sign(reg['uz_darcy_m_s']))

        if np.isfinite(reg['uz_darcy_m_s']) and mean_dh != 0 and signs_agree:
            K = darcy_law_K(reg['uz_darcy_m_s'], mean_dh, DZ_HEAD)
            k_intrinsic = intrinsic_permeability(K, mean_twater)
            log_k = np.log10(k_intrinsic)
            if verbose:
                print(f"Étape 5 - gradient dh/dz={mean_dh/DZ_HEAD:.3f} (dh moyen={mean_dh:.3f} m, "
                      f"{'stable' if gradient_stable else 'INSTABLE (change de signe dans la fenêtre)'}) "
                      f"-> K={K:.2e} m/s -> k={k_intrinsic:.2e} m² (log_k={log_k:.2f}, T_eau moy={mean_twater:.1f}°C)\n")
        else:
            K = k_intrinsic = log_k = np.nan
            if verbose:
                if np.isfinite(reg['uz_darcy_m_s']) and mean_dh != 0 and not signs_agree:
                    print(f"Étape 5 - DÉSACCORD DE SENS : dh moyen={mean_dh:.3f} m "
                          f"({'infiltration' if mean_dh > 0 else 'exfiltration'} attendue) mais uz thermique="
                          f"{reg['uz_darcy_mm_yr']:.0f} mm/an ({'infiltration' if reg['uz_darcy_mm_yr'] > 0 else 'exfiltration'}) "
                          "-> K non calculé (résultat non physique)\n")
                else:
                    print("Étape 5 - impossible (uz non défini ou gradient nul sur cette fenêtre)\n")

        reg['window'] = window_label
        reg['window_index'] = w_idx
        reg['t_start'] = t0 + pd.Timedelta(seconds=t_start)
        reg['t_end'] = t0 + pd.Timedelta(seconds=t_end)
        reg['mean_dh_m'] = mean_dh
        reg['gradient_stable'] = gradient_stable
        reg['signs_agree'] = signs_agree
        reg['K_m_s'] = K
        reg['k_intrinsic_m2'] = k_intrinsic
        reg['log_k'] = log_k
        rows.append({k: v for k, v in reg.items() if k not in ('names', 'z', 'log_amp', 'phase_unwrapped')})
        wdata[w_idx] = {'label': window_label, 'mask': mask, 'fits': fits, 'reg': reg}

    return pd.DataFrame(rows), wdata


def porosity_envelope(coherent_df, reliable_df, wdata):
    """
    Balaie n sur POROSITY_RANGE (même plage que le grid search) et renvoie
    l'enveloppe ((lam_min,lam_max), (log_k_min,log_k_max)) compatible avec
    cette incertitude, à partir des fenêtres cohérentes/fiables déjà calculées.
    Gamma_joint et v_thermal ne dépendent pas de n (cf. ÉTAPE 6) : seule leur
    conversion en lam/uz/K/k en dépend.
    """
    porosity_samples = np.linspace(*POROSITY_RANGE, 10)

    # La conductivité thermique de la fraction SOLIDE doit toujours être supérieure
    # à celle de l'eau (les minéraux conduisent mieux la chaleur que l'eau) : même le
    # k_bulk d'un milieu saturé ne peut pas être inférieur à kw (un mélange eau+solide
    # ne peut que conduire au moins aussi bien que l'eau pure). Une paire (fenêtre,
    # porosité) qui donne k_bulk ou ks < kw signale un Gamma_joint peu fiable pour
    # cette fenêtre : on l'exclut de l'enveloppe plutôt que d'accepter un résultat
    # non physique.
    lam_all = []
    n_rejected = 0
    for _, row in coherent_df.iterrows():
        gamma_joint = row['gamma_joint_1e-6'] * 1e-6
        for n in porosity_samples:
            cv_n = n * CW_WATER + (1 - n) * CS_SOLID
            k_bulk = gamma_joint * cv_n
            if k_bulk < KW_WATER:
                n_rejected += 1
                continue
            ks = solid_thermal_conductivity(k_bulk, n, KW_WATER)
            if ks < KW_WATER:
                n_rejected += 1
                continue
            lam_all.append(ks)
    if n_rejected:
        print(f"({n_rejected} combinaisons (fenêtre, porosité) écartées de l'enveloppe lam : "
              f"k_bulk ou ks < kw={KW_WATER} W/m/K, non physique)")

    log_k_all = []
    for _, row in reliable_df.iterrows():
        w_idx_row = int(row['window_index'])
        twater_window = df['TempMolo'].to_numpy()[wdata[w_idx_row]['mask']]
        mean_tw = float(np.nanmean(twater_window))
        for n in porosity_samples:
            cv_n = n * CW_WATER + (1 - n) * CS_SOLID
            uz_n = row['v_thermal_m_s'] * cv_n / CW_WATER
            K_n = darcy_law_K(uz_n, row['mean_dh_m'], DZ_HEAD)
            log_k_all.append(np.log10(intrinsic_permeability(K_n, mean_tw)))

    lam_all, log_k_all = np.array(lam_all), np.array(log_k_all)
    lam_range = (lam_all.min(), lam_all.max()) if len(lam_all) else (np.nan, np.nan)
    log_k_range = (log_k_all.min(), log_k_all.max()) if len(log_k_all) else (np.nan, np.nan)
    return lam_range, log_k_range


# =============================================================================
# ÉTAPE 6 - BALAYAGE DES TAILLES DE FENÊTRE
# On applique toute l'analyse pour chaque taille de fenêtre candidate, sans
# détail fenêtre par fenêtre (verbose=False), pour comparer objectivement.
# =============================================================================
print("=" * 80)
print("BALAYAGE DES TAILLES DE FENÊTRE")
print("=" * 80)

sweep_rows = []
analyses = {}  # on garde (summary, wdata) de chaque taille pour réutilisation
for wd in WINDOW_SIZES_TO_TEST:
    summ, wdata = analyze_window_size(wd, verbose=False)
    analyses[wd] = (summ, wdata)

    coh = summ[summ['coherence_rel_diff'] < COHERENCE_THRESHOLD]
    rel = coh[coh['gradient_stable'] & np.isfinite(coh['log_k'])]
    coherence_rate = len(coh) / len(summ) if len(summ) else np.nan
    agreement_rate = len(rel) / len(coh) if len(coh) else np.nan

    if not coh.empty:
        (lam_lo, lam_hi), (lk_lo, lk_hi) = porosity_envelope(coh, rel, wdata)
    else:
        (lam_lo, lam_hi), (lk_lo, lk_hi) = (np.nan, np.nan), (np.nan, np.nan)

    sweep_rows.append({
        'window_days': wd, 'n_windows': len(summ),
        'coherence_rate': coherence_rate, 'agreement_rate': agreement_rate,
        'lam_lo': lam_lo, 'lam_hi': lam_hi, 'log_k_lo': lk_lo, 'log_k_hi': lk_hi,
    })

sweep_df = pd.DataFrame(sweep_rows)
sweep_df.to_csv(RESULTS_DIR / "stallman_window_sweep.csv", sep=";", index=False)
with pd.option_context('display.float_format', '{:.2f}'.format):
    print(sweep_df.to_string(index=False))
print(f"\nTableau complet sauvegardé dans {RESULTS_DIR / 'stallman_window_sweep.csv'}\n")

# On choisit la taille de fenêtre qui maximise le taux d'accord uz/dh (agreement_rate),
# parmi celles offrant assez de fenêtres pour être statistiquement raisonnables.
candidates = sweep_df[sweep_df['n_windows'] >= MIN_WINDOWS_FOR_CANDIDATE].dropna(subset=['agreement_rate'])
if candidates.empty:
    candidates = sweep_df.dropna(subset=['agreement_rate'])
RECOMMENDED_WINDOW_DAYS = int(candidates.loc[candidates['agreement_rate'].idxmax(), 'window_days'])
print(f"-> Taille de fenêtre retenue : {RECOMMENDED_WINDOW_DAYS} jours "
      f"(meilleur taux d'accord uz/dh parmi les tailles avec >= {MIN_WINDOWS_FOR_CANDIDATE} fenêtres)\n")

# =============================================================================
# Détail complet (ÉTAPES 1 à 5) pour la taille de fenêtre retenue
# =============================================================================
print("=" * 80)
print(f"DÉTAIL COMPLET POUR LA TAILLE DE FENÊTRE RETENUE ({RECOMMENDED_WINDOW_DAYS} jours)")
print("=" * 80)
summary, window_data = analyze_window_size(RECOMMENDED_WINDOW_DAYS, verbose=True)
summary.to_csv(RESULTS_DIR / "stallman_diurnal.csv", sep=";", index=False)
print(f"Résumé (une ligne par fenêtre) sauvegardé dans {RESULTS_DIR / 'stallman_diurnal.csv'}\n")

coherent = summary[summary['coherence_rel_diff'] < COHERENCE_THRESHOLD]
print(f"{len(coherent)}/{len(summary)} fenêtres jugées cohérentes "
      f"(désaccord Gamma_ampl/Gamma_phase < {COHERENCE_THRESHOLD:.0%})")

reliable_k = pd.DataFrame()
if coherent.empty:
    best_idx = summary['coherence_rel_diff'].idxmin()
    print("Aucune fenêtre suffisamment cohérente : bornes non calculées. "
          "On illustre malgré tout la méthode sur la fenêtre la moins mauvaise.")
else:
    print(f"- vitesse de Darcy (uz) : {coherent['uz_darcy_mm_yr'].min():.0f} à {coherent['uz_darcy_mm_yr'].max():.0f} mm/an")
    best_idx = coherent['coherence_rel_diff'].idxmin()

    reliable_k = coherent[coherent['gradient_stable'] & np.isfinite(coherent['log_k'])]
    print(f"\n{len(reliable_k)}/{len(coherent)} fenêtres cohérentes ont EN PLUS un gradient stable "
          "et un accord de signe uz/dh (nécessaire pour que log_k soit interprétable) :")
    if reliable_k.empty:
        print("Aucune fenêtre ne remplit les deux conditions : pas de bornes log_k proposées.")
    else:
        (lam_lo, lam_hi), (lk_lo, lk_hi) = porosity_envelope(coherent, reliable_k, window_data)
        print(f"\nBornes suggérées pour le grid search (n balayée sur {POROSITY_RANGE}) :")
        print(f"- lam (fraction solide) : {lam_lo:.2f} à {lam_hi:.2f} W/m/K")
        print(f"- log_k (perméabilité intrinsèque) : {lk_lo:.2f} à {lk_hi:.2f} (k = {10**lk_lo:.1e} à {10**lk_hi:.1e} m²)")
        # On préfère illustrer les graphiques pédagogiques avec une fenêtre où TOUTE
        # la chaîne (étapes 1 à 5) est cohérente, plutôt que la meilleure fenêtre au
        # seul sens thermique (qui peut avoir un désaccord de signe non physique).
        best_idx = reliable_k['coherence_rel_diff'].idxmin()

best_w_idx = int(summary.loc[best_idx, 'window_index'])
best = window_data[best_w_idx]
print(f"\nFenêtre retenue pour les graphiques pédagogiques ci-dessous : {best['label']} "
      f"(la plus cohérente, désaccord={best['reg']['coherence_rel_diff']:.1%})")

# =============================================================================
# PÉRIODE DE CALIBRATION RECOMMANDÉE
# On cherche la plus longue série de fenêtres CONSÉCUTIVES à la fois cohérentes
# (Gamma_ampl ~ Gamma_phase) et en accord de sens (uz thermique ~ dh) : c'est la
# période la plus digne de confiance pour caler Ginette (2_run_real_case.py).
# =============================================================================
print("\n" + "=" * 80)
print("PÉRIODE DE CALIBRATION RECOMMANDÉE")
print("=" * 80)

reliable_idx = set(reliable_k['window_index']) if not reliable_k.empty else set()
all_idx = sorted(summary['window_index'])
runs = []
current_run = []
for idx in all_idx:
    if idx in reliable_idx:
        current_run.append(idx)
    else:
        if current_run:
            runs.append(current_run)
        current_run = []
if current_run:
    runs.append(current_run)

if not runs:
    print("Aucune période consécutive de fenêtres fiables trouvée avec la taille de "
          f"fenêtre retenue ({RECOMMENDED_WINDOW_DAYS}j). Vérifie deltaP ou baisse COHERENCE_THRESHOLD.")
else:
    best_run = max(runs, key=len)
    run_start = window_data[best_run[0]]['reg']['t_start']
    run_end = window_data[best_run[-1]]['reg']['t_end']
    run_days = (run_end - run_start).days
    mean_coherence = summary[summary['window_index'].isin(best_run)]['coherence_rel_diff'].mean()
    print(f"Période la plus longue où thermique et charge s'accordent en continu "
          f"({len(best_run)} fenêtre(s) de {RECOMMENDED_WINDOW_DAYS}j consécutive(s)) :")
    print(f"-> DU {run_start.date()} AU {run_end.date()} ({run_days} jours, "
          f"désaccord Gamma_ampl/Gamma_phase moyen={mean_coherence:.1%})")
    print(f"\nPour 0_boundary_conditions_real_case.py / 2_run_real_case.py, utiliser :")
    print(f"  date_simul_bg = pd.to_datetime(\"{run_start.strftime('%Y/%m/%d %H:%M:%S')}\")")
    print(f"  nb_day = {run_days}")

    # Recommandation sauvegardée pour config_lomos.py (FORCE_DATE_SIMUL_BG=False
    # la relit directement, plutôt que de la recopier à la main dans le config).
    period_file = RESULTS_DIR / "stallman_recommended_period.txt"
    with open(period_file, "w") as f:
        f.write(f"date_simul_bg {run_start.strftime('%Y/%m/%d %H:%M:%S')}\n")
        f.write(f"nb_day {run_days}\n")
    print(f"Période recommandée sauvegardée dans {period_file}")

# =============================================================================
# GRAPHIQUE 1 - Séries temporelles brutes + onde diurne ajustée, à chaque
# profondeur, sur la fenêtre retenue. Permet de VOIR l'amortissement et le
# retard avec la profondeur (au lieu de seulement le lire dans un tableau).
# =============================================================================
mask = best['mask']
t_window = df.index[mask]
w = 2 * np.pi / PERIOD

fig1, axes = plt.subplots(len(sensor_cols), 1, figsize=(10, 10), sharex=True)
for ax, name in zip(axes, sorted(DEPTHS, key=lambda n: DEPTHS[n])):
    raw = df[name].to_numpy()[mask]
    f = best['fits'][name]
    t_rel = time_s[mask]
    t_trend = (t_rel - t_rel.mean()) / 86400.0
    trend_poly = np.poly1d(np.polyfit(t_trend, raw, TREND_ORDER))  # tendance seule, pour l'affichage
    fitted = trend_poly(t_trend) + f['amplitude'] * np.cos(w * t_rel - f['phase'])

    ax.plot(t_window, raw, '.', color='0.6', markersize=2, label='mesures brutes')
    ax.plot(t_window, fitted, 'r-', linewidth=1.2, label='tendance + onde diurne ajustée')
    ax.set_ylabel(f"{name}\n(z={DEPTHS[name]*100:.0f} cm)\n°C")
    ax.legend(loc='upper right', fontsize=8)
axes[0].set_title(f"Onde diurne ajustée à chaque profondeur - fenêtre {best['label']}\n"
                   "(l'amplitude diminue et le signal se décale avec la profondeur)")
axes[-1].set_xlabel('Date')
plt.tight_layout()
plt.savefig(RESULTS_DIR / "stallman_fit_temporel.png", dpi=200)

# =============================================================================
# GRAPHIQUE 2 - Le cœur pédagogique de la méthode : ln(amplitude) et phase en
# fonction de la profondeur, avec les droites de régression dont les pentes
# donnent directement a et b (voir ÉTAPE 2 dans l'entête du script).
# =============================================================================
reg = best['reg']
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

z_fine = np.linspace(reg['z'].min(), reg['z'].max(), 50)
ax1.plot(reg['z'] * 100, reg['log_amp'], 'o', color='tab:blue')
ax1.plot(z_fine * 100, reg['slope_amp'] * z_fine + reg['intercept_amp'], '--', color='tab:blue')
ax1.set_xlabel('Profondeur z [cm]')
ax1.set_ylabel('ln(amplitude) [ln(°C)]')
ax1.set_title(f"Amortissement : pente = -a = {reg['slope_amp']:.3f} (R²={reg['r2_amp']:.3f})")

ax2.plot(reg['z'] * 100, reg['phase_unwrapped'], 'o', color='tab:red')
ax2.plot(z_fine * 100, reg['slope_phase'] * z_fine + reg['intercept_phase'], '--', color='tab:red')
ax2.set_xlabel('Profondeur z [cm]')
ax2.set_ylabel('Phase [rad]')
ax2.set_title(f"Déphasage : pente = b = {reg['slope_phase']:.3f} (R²={reg['r2_phase']:.3f})")

fig2.suptitle(f"Fenêtre {best['label']} — Gamma_ampl={reg['gamma_ampl_1e-6']:.2f}e-6, "
              f"Gamma_phase={reg['gamma_phase_1e-6']:.2f}e-6 m²/s, "
              f"uz={reg['uz_darcy_mm_yr']:.0f} mm/an")
plt.tight_layout()
plt.savefig(RESULTS_DIR / "stallman_regression_profondeur.png", dpi=200)

# =============================================================================
# GRAPHIQUE 3 - Sensibilité à la taille de fenêtre : coherence_rate et
# agreement_rate en fonction de WINDOW_DAYS (visualise le compromis discuté à
# l'ÉTAPE 6, et justifie objectivement le choix de RECOMMENDED_WINDOW_DAYS).
# =============================================================================
fig3, ax = plt.subplots(figsize=(8, 5))
ax.plot(sweep_df['window_days'], sweep_df['coherence_rate'] * 100, 'o-', label='% fenêtres cohérentes (Γampl≈Γphase)')
ax.plot(sweep_df['window_days'], sweep_df['agreement_rate'] * 100, 's-', label='% en plus en accord uz/dh')
ax.axvline(RECOMMENDED_WINDOW_DAYS, color='0.5', linestyle=':', label=f'retenue ({RECOMMENDED_WINDOW_DAYS}j)')
ax.set_xlabel('Taille de fenêtre [jours]')
ax.set_ylabel('%')
ax.set_title('Sensibilité à la taille de fenêtre')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(RESULTS_DIR / "stallman_window_sweep.png", dpi=200)
plt.show()
