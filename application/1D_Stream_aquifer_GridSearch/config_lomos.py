#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Config partagée entre 0_boundary_conditions_real_case.py, 1_define_grid_search.py,
2_run_real_case.py, 3_misfit.py et 4_plot_results.py - un seul endroit à changer
pour basculer de point d'observation, de config BC, de période, de bornes du
grid search ou de pipeline à exécuter (voir README.md).
@author: Agnès Rivière, Samuel Larance, Alexandrine Gesret
"""

from pathlib import Path

import numpy as np
import pandas as pd

# =============================================================================
# POINT D'OBSERVATION ET CONFIG BC
# =============================================================================
# Point d'observation (dossier dans OBS_point/) - détermine aussi le
# sous-dossier results/{POINT_NAME}/ où chaque script lit/écrit, pour que deux
# points différents ne s'écrasent jamais mutuellement.
POINT_NAME = "lomos230"

# Profondeur de chaque capteur sous le lit du ru [m], positif vers le bas -
# caractéristique de la station de terrain (dispositif LOMOS-mini : 10cm entre
# capteurs jusqu'à 40cm). Une autre station peut avoir un autre espacement/une
# autre profondeur max, tant que l'espacement reste UNIFORME (voir DZ_OBS
# ci-dessous - Ginette place ses points d'observation à pas fixe).
SENSOR_DEPTHS = {'TempMolo': 0.00, 'Temp1': 0.10, 'Temp2': 0.20, 'Temp3': 0.30, 'Temp4': 0.40}

# Capteur où deltaP (différence de charge rivière/piézomètre) est réellement
# mesuré - fixe, ne dépend pas de TOP_SENSOR/BOTTOM_SENSOR choisis ci-dessous.
DELTAP_SENSOR = 'Temp4'

# Capteurs utilisés comme CL haute/basse du domaine simulé, à choisir librement
# parmi SENSOR_DEPTHS. Deux presets déjà étudiés sur lomos230 (voir README.md) :
# - Config A (TOP_SENSOR='TempMolo', BOTTOM_SENSOR='Temp4') : domaine naturel
#   complet, TempMolo brut (eau libre) en CL haute.
# - Config C (TOP_SENSOR='Temp1', BOTTOM_SENSOR='Temp4') : évite TempMolo (eau
#   libre, pas la même grandeur physique que le sédiment), domaine réduit,
#   log_k mieux identifié (Pe plus grand sur un domaine plus court).
TOP_SENSOR = 'Temp1'
BOTTOM_SENSOR = 'Temp4'

# Domaine simulé relabellisé sur TOP_SENSOR (Z_TOP=0 toujours, PAS la
# profondeur absolue du capteur) : les équations de Ginette n'utilisent z que
# par différences (vérifié par relecture du code + algèbre à la main), donc ce
# décalage est physiquement neutre tant que les CL de charge sont dérivées
# indépendamment de la physique et pas juste translatées avec z (voir README.md).
# round() : la soustraction flottante (ex: 0.40-0.10) peut introduire un bruit
# binaire (-0.30000000000000004 au lieu de -0.30) qui fait déborder d'une
# maille les np.arange de generate_zone_parameters/initial_conditions - leurs
# marges de arrêt diffèrent légèrement, donc ce bruit les désynchronise.
DOMAIN_LENGTH = round(abs(SENSOR_DEPTHS[BOTTOM_SENSOR] - SENSOR_DEPTHS[TOP_SENSOR]), 6)
Z_TOP = 0.0
Z_BOTTOM = -DOMAIN_LENGTH

# Ginette place toujours exactement 4 points d'observation, équidistants de
# DZ_OBS sous TOP_SENSOR (E_parametre_bck.dat a 4 emplacements [cell1]-[cell4],
# figés côté Fortran - voir setup_ginette* dans Direct_model.py). DZ_OBS est
# donc dérivé de l'espacement RÉEL entre capteurs, qui doit être uniforme.
_depths_sorted = sorted(SENSOR_DEPTHS.values())
_spacings = np.diff(_depths_sorted)
if not np.allclose(_spacings, _spacings[0]):
    raise ValueError(f"SENSOR_DEPTHS doit être uniformément espacé (Ginette place "
                      f"ses 4 points d'observation à pas fixe) - espacements trouvés : {_spacings}")
DZ_OBS = float(_spacings[0])

# 3 des 4 points d'observation Ginette servent de capteurs de comparaison pour
# le misfit (3_misfit.py ne gère que misfit_1/2/3) - le 4e est ignoré (déjà le
# cas avant cette généralisation, ex: Config C). Si l'un tombe sur
# BOTTOM_SENSOR, la comparaison y est triviale (la maille compare à sa propre
# CL) : BOTTOM_SENSOR_TRIVIAL_INDEX indique lequel, pour que 3_misfit.py
# l'exclue de misfit_fair.
_depth_to_sensor = {round(d, 6): s for s, d in SENSOR_DEPTHS.items()}
COMPARISON_SENSORS = []
for _i in range(1, 4):
    _depth = round(SENSOR_DEPTHS[TOP_SENSOR] + DZ_OBS * _i, 6)
    if _depth not in _depth_to_sensor:
        raise ValueError(f"Aucun capteur à {_depth}m sous la surface (point d'observation "
                          f"{_i} attendu par Ginette) - vérifier SENSOR_DEPTHS/TOP_SENSOR/BOTTOM_SENSOR.")
    COMPARISON_SENSORS.append(_depth_to_sensor[_depth])

BOTTOM_SENSOR_TRIVIAL_INDEX = (COMPARISON_SENSORS.index(BOTTOM_SENSOR)
                                if BOTTOM_SENSOR in COMPARISON_SENSORS else None)

# =============================================================================
# PÉRIODE ET DISCRÉTISATION TEMPORELLE
# =============================================================================
# FORCE_DATE_SIMUL_BG=True : utilise DATE_SIMUL_BG/NB_DAY imposés ci-dessous.
# FORCE_DATE_SIMUL_BG=False : reprend la période recommandée par l'analyse
# Stallman (stallman_diffusivity.py - plus longue fenêtre continue où
# le sens d'écoulement thermique et hydraulique s'accordent). Si le fichier de
# résultat n'existe pas encore, le script est lancé automatiquement ici.
FORCE_DATE_SIMUL_BG = True

if FORCE_DATE_SIMUL_BG:
    DATE_SIMUL_BG = pd.to_datetime("2026/04/30 04:00:00")
    NB_DAY = 42
else:
    _app_dir = Path(__file__).resolve().parent
    _stallman_script = _app_dir / "stallman_diffusivity.py"
    _stallman_period_file = _app_dir / "results" / POINT_NAME / "stallman_recommended_period.txt"

    if not _stallman_period_file.exists():
        import subprocess
        import sys as _sys
        print(f"FORCE_DATE_SIMUL_BG=False : période non encore calculée, "
              f"lancement de {_stallman_script.name}...")
        subprocess.run([_sys.executable, str(_stallman_script)], check=True, cwd=str(_app_dir))

    if not _stallman_period_file.exists():
        raise FileNotFoundError(
            f"{_stallman_script.name} n'a trouvé aucune période fiable (voir sa sortie "
            "ci-dessus) - repasser à FORCE_DATE_SIMUL_BG=True et fixer DATE_SIMUL_BG/NB_DAY "
            "à la main."
        )

    _stallman_period = dict(
        line.split(maxsplit=1) for line in _stallman_period_file.read_text().splitlines()
    )
    DATE_SIMUL_BG = pd.to_datetime(_stallman_period["date_simul_bg"])
    NB_DAY = int(_stallman_period["nb_day"])

DT = 900  # pas de temps [s], doit correspondre au pas des observations

# =============================================================================
# RESSOURCES CALCUL
# =============================================================================
MAX_WORKERS = 6  # poste local, pas un cluster - laisser de la marge

# Sauvegarde de sim_velocity_{ID}.txt (vitesse Ginette réelle, ~9.5 Mo/simulation,
# ~5.4 Go pour un grid search à 576 points) : 3_misfit.py ne s'en sert que pour la
# meilleure simulation (voir ginette_velocity_for_row), avec repli automatique sur
# la vitesse Darcy analytique si absent. False par défaut pour éviter la charge
# disque/I/O d'un run complet (a contribué à un plantage machine le 2026-07-24) -
# passer à True seulement si on veut inspecter les profils de vitesse de TOUTES
# les simulations, pas juste la meilleure.
SAVE_VELOCITY_PROFILES = False

# =============================================================================
# PIPELINE : quelles étapes lancer depuis 2_run_real_case.py
# =============================================================================
# RUN_GRID_SEARCH=False : ne relance pas les simulations Ginette, réutilise
# celles déjà présentes dans results/{POINT_NAME}/ (utile pour ne refaire que
# le misfit/les plots après un calage déjà terminé).
RUN_GRID_SEARCH = True
# RUN_MISFIT=True : enchaîne automatiquement sur 3_misfit.py après le grid search.
RUN_MISFIT = True
# RUN_PLOTS=True : enchaîne automatiquement sur 4_plot_results.py après le misfit.
RUN_PLOTS = True

# =============================================================================
# BORNES DU GRID SEARCH (1_define_grid_search.py)
# =============================================================================
# Nombre de zones géologiques de la colonne : 1 = homogène, 2 = deux couches
# séparées par une interface à ALT_THK (zone 2 = couche du haut). NB_ZONE=1
# retenu suite au diagnostic TempMolo - historique dans la mémoire projet.
NB_ZONE = 1
ALT_THK = -0.05  # profondeur de l'interface [m] si NB_ZONE=2 (sans effet sinon)

# Nombre de valeurs testées par paramètre : repli (np.linspace) pour tout
# paramètre absent de PARAM_STEP ci-dessous.
N = 8

# Paramètres réellement testés dans le grid search : chaque nom doit
# correspondre à une variable [min, max] définie ci-dessous.
Name_parameters = ["log_k", "lam", "n", "c"]

# Bornes [min, max] - large par défaut (>= 4 ordres de grandeur sur log_k)
# pour ne pas présupposer la valeur calée et pouvoir estimer l'incertitude.
log_k = [-18, -11]  # log10(perméabilité intrinsèque k [m2])
lam = [1, 6]        # conductivité thermique de la fraction solide [W/m/K]
n = [5e-2, 0.6]     # porosité [-]
# c = capacité calorifique SPÉCIFIQUE (pas volumique) du solide [J/kg/K],
# combinée à rhosi (E_p_therm.dat, fixe) via CASE('ZHZ'). rhosi a été baissé à
# 1180 kg/m3 (sédiment riche en matière organique, pas du quartz pur) pour que
# c=2000, mélangé à l'eau selon la porosité (Cv = n*Cw_vol + (1-n)*rhos*c),
# donne une capacité volumique du milieu SATURÉ cohérente avec la littérature
# streambed (~2-3.5 MJ/m3/K, Stallman/Lapham/Constantz) - voir E_p_therm_bck.dat.
c = [2000]          # capacité thermique spécifique de la fraction solide [J/kg/K]

# Pas fixe (résolution) par paramètre : {nom: (pas, décimales d'arrondi)}.
# La grille est alors np.arange(min, max+pas, pas) arrondie à `décimales`. Un
# paramètre de Name_parameters absent d'ici retombe sur N/np.linspace.
# Pas resserrés à un moment (log_k=0.5, lam=0.1) sans revérifier la taille
# totale -> 9180 combinaisons, très au-dessus de la limite de 500 sur ce
# poste (voir mémoire feedback_grid_search_cpu_limit). Remis à une résolution
# plus grossière (proche des pas d'origine) : 8x6x7x1 = 336 combinaisons.
PARAM_STEP = {
    "log_k": (1, 1),
    "lam": (1, 1),
    "n": (0.1, 3),
}