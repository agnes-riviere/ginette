#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Config partagée entre 0_boundary_conditions_real_case.py, 1_define_grid_search.py,
2_run_real_case.py, 3_misfit.py et 4_plot_results.py - un seul endroit à changer
pour basculer de point d'observation, de config BC, de période, de bornes du
grid search ou de pipeline à exécuter (voir README.md).
@author: Agnès Rivière, Samuel Larance, Alexandrine Gesret
"""

import os
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
# Surcharge par l'environnement (balayage d'offset sur un autre point sans
# éditer la config - 5_offset_sweep.py).
POINT_NAME = os.environ.get("GINETTE_POINT_NAME", POINT_NAME)

# Profondeur de chaque capteur sous le lit du ru [m], positif vers le bas -
# caractéristique de la station de terrain (dispositif LOMOS-mini : 10cm entre
# capteurs jusqu'à 40cm). Une autre station peut avoir un autre espacement/une
# autre profondeur max, tant que l'espacement reste UNIFORME (voir DZ_OBS
# ci-dessous - Ginette place ses points d'observation à pas fixe).
SENSOR_DEPTHS = {'TempMolo': 0.00, 'Temp1': 0.10, 'Temp2': 0.20, 'Temp3': 0.30, 'Temp4': 0.40}

# Capteur où deltaP (différence de charge rivière/piézomètre) est réellement
# mesuré - fixe, ne dépend pas de TOP_SENSOR/BOTTOM_SENSOR choisis ci-dessous.
DELTAP_SENSOR = 'Temp4'
# Inclinaison du tube de pression par rapport à la verticale [deg], par point.
# Le capteur étant différentiel (deux prises en eau sur le même diaphragme),
# l'angle ne fausse pas dp ; il rapproche seulement la crépine de la surface :
# profondeur verticale = longueur x cos(angle). lomos231 : tube posé à 16°
# (info terrain 2026-09-15), crépine à 0.40 x cos(16°) = 0.385 m -> gradient
# +3.9 %, soit Δlog_k ≈ -0.02, invisible au pas de grille. Points absents : 0°.
DELTAP_TUBE_ANGLE_DEG = {'lomos231': 16.0}
DELTAP_DEPTH = SENSOR_DEPTHS[DELTAP_SENSOR] * np.cos(np.radians(DELTAP_TUBE_ANGLE_DEG.get(POINT_NAME, 0.0)))

# Convention de dp (H_corrige_24h) : dp = h_riviere - h_nappe, donc dp > 0 =
# infiltration, dp < 0 = exfiltration. Le signe est celui des données et
# n'est jamais modifié ici ; seul DELTAP_OFFSET (décalage de zéro du capteur,
# ci-dessous) peut le faire basculer.

# Offset de zéro du capteur de pression différentielle [m], ajouté à deltaP
# APRÈS conversion en mètres et avant la réduction au sous-domaine (x0.75 en
# config C). 0 = capteur supposé juste. Sur lomos231, dp ne fait que +2 à
# +3 cm à l'étiage : un décalage de zéro de cet ordre suffit à inverser le
# sens réel de l'écoulement. Valeur pour lomos231 CALÉE (pas mesurée, aucune
# mesure manuelle rivière/piézo disponible) : balayage -2 à -7 cm sur la
# période 12/03->30/05/2026, misfit complet minimal à -4 cm, misfit débiaisé
# minimal à -5 cm, log_k=-11.5 dans les deux cas -> -4.5 +/- 1 cm. Sans
# offset (0), le calage part en bord de grille (log_k<=-15, lam=1, n=0.65)
# car les températures indiquent une exfiltration que dp>0 interdit.
# Par point : les autres points (ex: lomos230, tube vertical) restent à 0 -
# le même balayage sur lomos230 (3 périodes) ne réclame aucun offset et
# rejette celui du 231 (misfit dégradé dès -4 cm).
DELTAP_OFFSET_CM = {'lomos231': -3.8}
DELTAP_OFFSET = DELTAP_OFFSET_CM.get(POINT_NAME, 0.0) / 100.0
# Surcharge temporaire par 5_offset_sweep.py (balayage de l'offset, lance le
# pipeline une fois par valeur) - ne rien changer ici pour un balayage.
DELTAP_OFFSET = float(os.environ.get("GINETTE_DELTAP_OFFSET", DELTAP_OFFSET))

# Test de sensibilité : Δh CONSTANT imposé sur toute la période [cm], à la place
# de la série dp mesurée (déjà ramené au sous-domaine : c'est h_top - h_bottom
# tel que Ginette le voit). None = série mesurée (défaut). Sert à voir si la
# variabilité temporelle de dp apporte quelque chose au calage.
DELTAP_CONST_CM = None
_env_const = os.environ.get("GINETTE_DELTAP_CONST_CM")
if _env_const not in (None, "", "None"):
    DELTAP_CONST_CM = float(_env_const)

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

# Capteurs de CALAGE = ceux strictement entre TOP_SENSOR et BOTTOM_SENSOR :
# les deux capteurs de CL ne sont pas comparés (la maille de CL reproduit sa
# propre CL, la comparaison y serait triviale). Config C (Temp1/Temp4) ->
# Temp2 et Temp3 ; Config A (TempMolo/Temp4) -> Temp1, Temp2, Temp3.
# Ginette place ses points d'observation à DZ_OBS, 2*DZ_OBS, 3*DZ_OBS sous la
# CL haute ([cell1]..[cell3] de E_parametre_bck.dat, lus par run_direct_model) :
# le capteur à la profondeur p correspond au point n° (p - p_top)/DZ_OBS.
# Les noms de capteurs PHYSIQUES sont conservés de bout en bout
# (observed_data.txt, sim_temp_ID.txt, colonnes misfit_2/misfit_3... de
# results.txt) - jamais renumérotés.
CALIB_SENSORS = {}   # nom du capteur -> n° du point d'observation Ginette (1..3)
for _s, _p in sorted(SENSOR_DEPTHS.items(), key=lambda kv: kv[1]):
    if SENSOR_DEPTHS[TOP_SENSOR] < _p < SENSOR_DEPTHS[BOTTOM_SENSOR]:
        _k = (_p - SENSOR_DEPTHS[TOP_SENSOR]) / DZ_OBS
        if abs(_k - round(_k)) > 1e-6 or not 1 <= round(_k) <= 3:
            raise ValueError(f"Capteur {_s} à {_p} m : pas de point d'observation Ginette à "
                             f"cette profondeur (points à 1, 2, 3 x DZ_OBS={DZ_OBS} sous {TOP_SENSOR}).")
        CALIB_SENSORS[_s] = int(round(_k))
if not CALIB_SENSORS:
    raise ValueError(f"Aucun capteur entre {TOP_SENSOR} et {BOTTOM_SENSOR} : rien à caler.")


def sensor_index(name):
    """Indice utilisé dans les colonnes de results.txt (misfit_2, mae_2, ...) :
    le numéro du capteur PHYSIQUE ('Temp2' -> '2'), pas un rang."""
    return name.replace("Temp", "")


# =============================================================================
# PÉRIODE ET DISCRÉTISATION TEMPORELLE
# =============================================================================
# DATE_SIMUL_BG/NB_DAY = fenêtre COMPARÉE aux observations (misfit). La
# simulation elle-même démarre SPIN_UP_RUN_DAYS jours plus tôt (voir la section
# MISFIT plus bas : RUN_DATE_BG/RUN_NB_DAY).
# FORCE_DATE_SIMUL_BG=True : utilise DATE_SIMUL_BG/NB_DAY imposés ci-dessous.
# FORCE_DATE_SIMUL_BG=False : reprend la période recommandée par l'analyse
# Stallman (stallman_diffusivity.py - plus longue fenêtre continue où
# le sens d'écoulement thermique et hydraulique s'accordent). Si le fichier de
# résultat n'existe pas encore, le script est lancé automatiquement ici.
FORCE_DATE_SIMUL_BG = True

# Période PAR POINT (début de la fenêtre comparée, nb de jours) : évite de
# simuler un point avec la période calée pour l'autre quand on change
# POINT_NAME. Un point absent du dict retombe sur DEFAULT_PERIOD.
PERIODS = {
    'lomos231': ("2026/03/12 00:00:00", 79),   # étiage + réchauffement, offset dp calé dessus
    'lomos230': ("2026/04/30 04:00:00", 69),   # 30/04 -> 08/07 (fin des données)
}
DEFAULT_PERIOD = ("2026/04/30 04:00:00", 42)

if FORCE_DATE_SIMUL_BG:
    _date_bg, _nb_day = PERIODS.get(POINT_NAME, DEFAULT_PERIOD)
    DATE_SIMUL_BG = pd.to_datetime(_date_bg)
    NB_DAY = int(_nb_day)
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
# Nombre fixe (ex: 6), ou "auto" = coeurs disponibles - 2. La surcharge
# d'environnement permet de limiter les workers sans modifier le dépôt.
_max_workers_env = os.environ.get("GINETTE_MAX_WORKERS", "auto").strip().lower()
try:
    MAX_WORKERS = max(1, int(_max_workers_env))
except ValueError:
    MAX_WORKERS = "auto"
# disponibles sur la machine qui lance le run, moins 2 (voir
# _n_worker_processes() dans 2_run_real_case.py) - pratique si on change
# souvent de poste/cluster sans repasser ici à chaque fois.

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
# Surcharge par 5_offset_sweep.py (pas de figures à chaque valeur d'offset).
RUN_PLOTS = bool(int(os.environ.get("GINETTE_RUN_PLOTS", int(RUN_PLOTS))))
# Les fenêtres graphiques sont désactivées par défaut pour que le pipeline soit
# sûr en batch, sur cluster et dans les workers multiprocessing. Activer avec
# GINETTE_SHOW_PLOTS=1 pour un lancement interactif local.
SHOW_PLOTS = bool(int(os.environ.get("GINETTE_SHOW_PLOTS", "0")))
# PLOT_BC=True : 0_boundary_conditions_real_case.py affiche les figures de
# contrôle (profils initiaux, CL de température/charge vs données brutes).
# False pour un enchaînement sans interruption (batch, cluster, balayage).
PLOT_BC = SHOW_PLOTS
PLOT_BC = bool(int(os.environ.get("GINETTE_PLOT_BC", int(PLOT_BC))))

# =============================================================================
# BALAYAGE DE L'OFFSET DU CAPTEUR DE PRESSION (5_offset_sweep.py)
# =============================================================================
# Valeurs de DELTAP_OFFSET testées [cm] : pour chacune, le pipeline complet
# (CL -> grille -> simulations -> misfit) est relancé et le meilleur misfit de
# la grille est relevé. Résultats dans results/{POINT_NAME}/offset_sweep/.
OFFSET_SWEEP_CM = list(range(-8, 3))  # 231 : encadre les -4.5 cm calés et le zéro

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

# Bornes [min, max] par point - large par défaut (>= 4 ordres de grandeur sur
# log_k) pour ne pas présupposer la valeur calée et pouvoir estimer
# l'incertitude. lomos230 : sable gréseux. lomos231 : hypothèse d'un lit
# granitique fracturé (bornes plus resserrées).
LOG_K_BY_POINT = {
    "lomos230": [-13, -11],  # sable gréseux : log10(perméabilité intrinsèque k [m2])
    "lomos231": [-15, -11],  # granite fracturé
}
LAM_BY_POINT = {
    "lomos230": [2, 6],      # sable gréseux : conductivité de la fraction solide [W/m/K]
    "lomos231": [2.5, 4.0],  # granite fracturé
}
N_BY_POINT = {
    "lomos230": [0.15, 0.65],  # sable gréseux : porosité totale plausible [-]
    "lomos231": [0.20, 0.40],  # granite fracturé
}
log_k = LOG_K_BY_POINT.get(POINT_NAME, LOG_K_BY_POINT["lomos230"])
lam = LAM_BY_POINT.get(POINT_NAME, LAM_BY_POINT["lomos230"])
n = N_BY_POINT.get(POINT_NAME, N_BY_POINT["lomos230"])
# c = capacité calorifique SPÉCIFIQUE (pas volumique) du solide [J/kg/K],
# combinée à rhosi (E_p_therm.dat, fixe) via CASE('ZHZ'). rhosi a été baissé à
# 1180 kg/m3 (sédiment riche en matière organique, pas du quartz pur) pour que
# c=2000, mélangé à l'eau selon la porosité (Cv = n*Cw_vol + (1-n)*rhos*c),
# donne une capacité volumique du milieu SATURÉ cohérente avec la littérature
# streambed (~2-3.5 MJ/m3/K, Stallman/Lapham/Constantz) - voir E_p_therm_bck.dat.

# Propriétés du solide, définies au même endroit que les bornes de calibration.
# lomos230 : sable gréseux riche en matière organique -> rhosi=1180 et c=2000
# J/kg/K déjà retenus lors du calage précédent (donnent une capacité volumique
# du milieu saturé cohérente avec la littérature streambed, ~2-3.5 MJ/m3/K,
# Stallman/Lapham/Constantz - voir E_p_therm_bck.dat). lomos231 : hypothèse
# d'un lit granitique fracturé.
RHO_SOLID_BY_POINT = {
    "lomos230": 1180.0,  # sable gréseux [kg/m3]
    "lomos231": 2650.0,  # granite [kg/m3]
}
C_SOLID_BY_POINT = {
    "lomos230": 2000.0,  # sable gréseux [J/kg/K]
    "lomos231": 800.0,   # granite [J/kg/K]
}
RHO_SOLID = RHO_SOLID_BY_POINT.get(POINT_NAME, RHO_SOLID_BY_POINT["lomos230"])
C_SOLID = C_SOLID_BY_POINT.get(POINT_NAME, C_SOLID_BY_POINT["lomos230"])
c = [C_SOLID]        # capacité thermique spécifique de la fraction solide [J/kg/K] (valeur fixe, pas calée)
# Pas fixe (résolution) par paramètre : {nom: (pas, décimales d'arrondi)}.
# La grille est alors np.arange(min, max+pas, pas) arrondie à `décimales`. Un
# paramètre de Name_parameters absent d'ici retombe sur N/np.linspace.
# Pas resserrés à un moment (log_k=0.5, lam=0.1) sans revérifier la taille
# totale -> 9180 combinaisons, très au-dessus de la limite de 500 sur ce
# poste (voir mémoire feedback_grid_search_cpu_limit). Remis à une résolution
# plus grossière (proche des pas d'origine) : 8x6x7x1 = 336 combinaisons.
PARAM_STEP = {
    "log_k": (0.25, 2),
    "lam": (1, 1),
    "n": (0.1, 3),
}

# =============================================================================
# MISFIT (3_misfit.py)
# =============================================================================
# Spin-up [jours] : le temps que le modèle "oublie" l'état initial (profil
# interpolé sur un seul instant, voir 0_boundary_conditions_real_case.py). Il
# est AJOUTÉ AVANT DATE_SIMUL_BG : la simulation démarre SPIN_UP_RUN_DAYS jours
# plus tôt (RUN_DATE_BG) et dure NB_DAY + SPIN_UP_RUN_DAYS (RUN_NB_DAY), et
# 3_misfit.py exclut ces premiers jours de la comparaison. La fenêtre
# DATE_SIMUL_BG -> DATE_SIMUL_BG + NB_DAY est donc exactement celle comparée
# aux observations, identique pour toutes les simulations de la grille. Les
# observations doivent exister SPIN_UP_RUN_DAYS jours avant DATE_SIMUL_BG
# (vérifié par 0_boundary_conditions_real_case.py).
#   entier : SPIN_UP_RUN_DAYS = SPIN_UP_DAYS.
#   "auto" : max(SPIN_UP_MIN_DAYS, SPIN_UP_TAU_MULTIPLIER x tau), tau = L^2/Gamma
#            le temps de diffusion thermique du matériau LE MOINS diffusif de la
#            grille (lam min, n max, c min) - un spin-up suffisant pour lui l'est
#            pour tous les autres. 3 tau = ~95% de l'état initial oublié
#            (relaxation exponentielle, exp(-3)), 5 tau = ~99%. Sur la grille
#            230 (L=0.3 m) : tau=5 j pour lam=1/n=0.65 -> 3 tau = 16 j, 5 tau = 26 j.
SPIN_UP_DAYS = 15
SPIN_UP_TAU_MULTIPLIER = 3
SPIN_UP_MIN_DAYS = 0

# Densité du solide [kg/m3], fixe (E_p_therm_bck.dat, rhosi) - voir le
# commentaire de `c` plus haut. Utilisée pour Cv dans 3_misfit.py et ici pour tau.
def _spin_up_run_days():
    if SPIN_UP_DAYS != "auto":
        return max(1, int(SPIN_UP_DAYS))
    import sys as _sys
    _src_py = str(Path(__file__).resolve().parents[2] / "src" / "src_python")
    if _src_py not in _sys.path:
        _sys.path.insert(0, _src_py)
    from Analytical_validation import bulk_thermal_conductivity, volumetric_heat_capacity
    lam_min, n_max, c_min = min(lam), max(n), min(c)
    gamma = (bulk_thermal_conductivity(lam_min, n_max)
             / volumetric_heat_capacity(n_max, RHO_SOLID, cs_specific=c_min))
    tau_days = DOMAIN_LENGTH**2 / gamma / 86400
    return max(1, int(SPIN_UP_MIN_DAYS), int(np.ceil(SPIN_UP_TAU_MULTIPLIER * tau_days)))


SPIN_UP_RUN_DAYS = _spin_up_run_days()
# Période effectivement simulée (0_boundary_conditions_real_case.py, 2_run_real_case.py).
# Surcharge par l'environnement (5_offset_sweep.py sur une autre période que
# celle de la config, pour valider l'offset sur des fenêtres indépendantes).
if os.environ.get("GINETTE_DATE_SIMUL_BG"):
    DATE_SIMUL_BG = pd.to_datetime(os.environ["GINETTE_DATE_SIMUL_BG"])
if os.environ.get("GINETTE_NB_DAY"):
    NB_DAY = int(os.environ["GINETTE_NB_DAY"])

RUN_DATE_BG = DATE_SIMUL_BG - pd.Timedelta(days=SPIN_UP_RUN_DAYS)

# Signature des CL que 0_boundary_conditions_real_case.py écrit dans
# GINETTE_SENSI/ (fichier bc_manifest.txt). 2_run_real_case.py la compare à la
# config courante avant de lancer : GINETTE_SENSI/ est un dossier de travail
# UNIQUE partagé entre les points, et sans ce contrôle un changement de
# POINT_NAME (ou de période/offset) sans repasser par l'étape 0 fait tourner
# un point avec les CL de l'autre, silencieusement.
def bc_manifest():
    return {
        "point": POINT_NAME,
        "run_date_bg": str(RUN_DATE_BG),
        "run_nb_day": str(RUN_NB_DAY),
        "top_sensor": TOP_SENSOR,
        "bottom_sensor": BOTTOM_SENSOR,
        "deltap_offset_m": f"{DELTAP_OFFSET:.5f}",
        "deltap_const_cm": str(DELTAP_CONST_CM),
        "deltap_depth_m": f"{DELTAP_DEPTH:.4f}",
    }
RUN_NB_DAY = NB_DAY + SPIN_UP_RUN_DAYS

# Erreur totale estimée sur la température (mesure + representativite du
# modele 1D), utilisée pour pondérer misfit_L2/misfit_L1 - pas juste la
# precision instrumentale (0.25 degC).
ERR_MISFIT = 5