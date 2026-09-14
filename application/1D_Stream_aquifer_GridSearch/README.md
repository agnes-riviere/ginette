---
noteId: "66ea7d50876b11f18009219d80a29152"
tags: []

---

# Calage de Ginette sur données réelles LOMOS-mini (lomos230)

Ce dossier calibre Ginette sur des mesures de terrain (dispositif LOMOS-mini, Cucchi et al.
2018 : 4 sondes de température à 10/20/30/40cm de profondeur dans le lit d'une rivière, un
capteur de température en eau libre TempMolo, et un capteur de différence de charge deltaP)
pour estimer les échanges nappe-rivière (perméabilité, porosité, conductivité thermique du
streambed).

## Pipeline

1. `0_boundary_conditions_real_case.py` — lit le CSV d'observation, construit les conditions
   initiales et aux limites, écrit les fichiers d'entrée Ginette (`GINETTE_SENSI/`).
2. `1_define_grid_search.py` — génère `results/grid_search.csv`, la grille de paramètres à
   tester (log_k, lam, n, c).
3. `2_run_real_case.py` — lance une simulation Ginette par ligne de la grille (parallélisé,
   `MAX_WORKERS=6` pour ne pas saturer un poste de travail local).
4. `3_misfit.py` — compare chaque simulation aux observations, calcule les métriques d'ajustement
   (misfit L2, MAE, KGE, pbias) et un diagnostic d'identifiabilité (nombre de Péclet, largeur du
   plateau de misfit).
5. `4_plot_results.py` — figures de synthèse.

Tous les paramètres utilisateur (point d'observation, config BC A/C, période de simulation,
`MAX_WORKERS`, `SAVE_VELOCITY_PROFILES`, bornes et pas du grid search) sont centralisés dans
`config_lomos.py` — un seul fichier à éditer, partagé par les 5 scripts du pipeline et
`stallman_diffusivity.py`.

La période de simulation (`DATE_SIMUL_BG`/`NB_DAY`) peut être imposée à la main
(`FORCE_DATE_SIMUL_BG=True`, valeurs ci-dessous) ou reprise automatiquement de la
recommandation de `stallman_diffusivity.py` (`FORCE_DATE_SIMUL_BG=False` : relit
`results/stallman_recommended_period.txt` ; si ce fichier n'existe pas encore, `config_lomos.py`
lance lui-même `stallman_diffusivity.py` pour le générer).

`stallman_diffusivity.py` lit `POINT_NAME` depuis `config_lomos.py` comme les autres scripts et
marche donc sur n'importe quel point au format LOMOS-mini (lomos230, lomos231, ...) ; seules les
profondeurs de capteurs (10/20/30/40cm) restent codées en dur, propres au dispositif LOMOS-mini.

`config_lomos.py` contrôle aussi l'enchaînement à partir de `2_run_real_case.py` :
- `RUN_GRID_SEARCH` : relance les simulations Ginette (False = réutilise les `sim_temp_*.txt`
  déjà présents, utile pour refaire juste le misfit/les plots).
- `RUN_MISFIT` : enchaîne sur `3_misfit.py` après le grid search.
- `RUN_PLOTS` : enchaîne sur `4_plot_results.py` après le misfit.

Chaque script ne déclenche que l'étape suivante (2→3→4), donc lancer `3_misfit.py` seul avec
`RUN_PLOTS=True` régénère aussi les figures, sans repasser par `2_run_real_case.py`.

## Le problème de la condition limite haute : TempMolo n'est pas la température du sédiment

TempMolo est mesuré 1-2cm au-dessus du lit, en eau libre — ce n'est pas la même grandeur
physique que la température du milieu poreux à z=0 (mélange convectif de l'eau libre vs
conduction/advection dans le sédiment). L'utiliser tel quel comme condition de Dirichlet en
surface revient à imposer au modèle un forçage dont il ne peut pas reproduire correctement la
physique, quels que soient les paramètres testés.

Le choix des capteurs de CL vit dans `config_lomos.py` (`TOP_SENSOR`/`BOTTOM_SENSOR`, librement
choisis parmi `SENSOR_DEPTHS`) : domaine, points de comparaison, rescaling de deltaP en découlent
tous automatiquement (voir la section Pipeline plus haut). Deux presets déjà étudiés sur lomos230 :

- **Config A** (`TOP_SENSOR='TempMolo'`, `BOTTOM_SENSOR='Temp4'`) : TempMolo brut en
  CL haute, Temp4 en CL basse, domaine 0 à -0.4m, calage sur Temp1/Temp2/Temp3.
- **Config C** (`TOP_SENSOR='Temp1'`, `BOTTOM_SENSOR='Temp4'`, par défaut, recommandée) : Temp1 ET Temp4, mesures réelles
  sans aucune correction, en CL haute/basse, domaine 0 à -0.3m, calage sur Temp2/Temp3 seulement.
  deltaP est ramené au sous-domaine (x0.75) en supposant un gradient hydraulique homogène sur la
  colonne.

Le nombre de points d'observation Ginette reste figé à 4, équidistants de `DZ_OBS` sous
`TOP_SENSOR` (format `E_parametre_bck.dat`, partagé avec d'autres apps) : `DZ_OBS` est dérivé de
`SENSOR_DEPTHS`, qui doit donc rester uniformément espacé (ex: 10cm entre capteurs pour
LOMOS-mini) - seule la profondeur/l'espacement peuvent varier d'une station à l'autre, pas le
nombre de capteurs.

Une troisième option (TempMolo corrigé par extrapolation Stallman depuis la paire fiable
Temp1-Temp2) a été testée puis abandonnée : la correction ne touche que l'harmonique 24h,
laissant la tendance sur plusieurs jours (dominante sur 42 jours de calage) inchangée. Étendre
la correction à tout le spectre par FFT est mal posé — le gain de Stallman croît avec la racine
de la fréquence et amplifie le bruit haute fréquence de façon exponentielle, produisant des
températures non physiques même avec un filtre passe-bas.

## Config A et Config C ne sont pas juste "différentes" : elles n'identifient pas la même chose

Le nombre de Péclet thermique (Pe = flux advectif / flux conductif, Kurylyk et al. 2017 eq.12)
dit a priori quels paramètres sont identifiables à partir de la seule température (Cucchi,
Flipo, Rivière, Rubin 2021, qui utilise Ginette sur ce même type de site) :
- Pe < 0.5 (conductif) : seule la diffusivité thermique effective est identifiable, la
  perméabilité seulement bornée.
- 0.5 < Pe < 5 (transition) : perméabilité ET diffusivité identifiables.
- Pe > 5 (advectif) : perméabilité identifiable, diffusivité non identifiable.

Sur lomos230 (42 jours, 2026-04-30 à 2026-06-11) :

| | log_k best-fit | \|Pe\| | plateau top-5% du misfit sur log_k | interprétation |
|---|---|---|---|---|
| Config A | -16 | ≈2×10⁻⁴ | [-18, -14] (4 ordres de grandeur) | **log_k non identifiable** — "-16" est un point pris au hasard sur un plateau plat |
| Config C | -12 | ≈0.63 | {-12} (une seule valeur de grille) | log_k bien identifié |

Config C fit aussi mieux sur les capteurs communs, en comparaison honnête (le `kge_3` de C
compare la maille de bord à sa propre condition limite, donc trivial, à exclure) : KGE
Temp2/Temp3 = 0.94/0.96 pour C contre 0.84/0.93 pour les mêmes profondeurs en A. Et le log_k=-12
de Config C est cohérent avec l'estimation Stallman indépendante (analyse spectrale pure, sans
Ginette) : -13.4 à -12.1 (`stallman_diffusivity.py`).

*Chiffres recalculés après les corrections du 2026-07-24 (voir section suivante) ; à
reconfirmer si tu relances un grid search complet, notamment `log_k` best-fit qui dépend
indirectement de Cv via le spin-up automatique.*

**Config C est recommandée par défaut sur ce site.** Config A reste disponible mais son log_k
doit être lu comme une borne, pas une estimation. lam et n restent mal contraints
individuellement même en Config C (équifinalité) — seule leur combinaison (diffusivité
effective) semble contrainte par les données ; éviter de sur-interpréter leurs valeurs
ponctuelles.

## Grille de paramètres

Volontairement large et non-informative même là où on a une idée a priori des paramètres
(log_k sur 7 ordres de grandeur : -18 à -11), pour pouvoir caractériser l'incertitude/le
plateau de misfit plutôt que de la restreindre a priori autour d'une estimation.

`c` (capacité calorifique **spécifique** du solide, J/kg/K, fixe à 2000 - PAS volumique malgré
ce que suggère son nom court) est combinée à `rhosi` (densité du solide, `E_p_therm.dat`, fixe
à 1180 kg/m3) via `CASE('ZHZ')` (`src/ginette_V2.f90`). 1180 kg/m3 (baissé depuis 2650/quartz
pur le 2026-07-24) reflète un sédiment riche en matière organique - cohérent avec c=2000, qui
est élevé pour un solide purement minéral. Ensemble, mélangés à l'eau selon la porosité
(Cv = n\*Cw_vol + (1-n)\*rhosi\*c), ils donnent une capacité volumique du milieu saturé de
2.45 à 3.46 MJ/m3/K sur n=[0.05,0.6], cohérente avec la littérature streambed (~2-3.5 MJ/m3/K,
Stallman/Lapham/Constantz).

## Convention de signe du flux (`darcy_flux_m_s`, `peclet` dans `3_misfit.py`)

Positif = vers le HAUT (exfiltration), négatif = vers le BAS (infiltration) - convention de
Ginette elle-même (z positif vers le haut, `src/ginette_V2.f90` ~l.3712 :
`v = -K/mu*(dP/dz+rho*g)`), pas celle de Stallman (z=profondeur, positif=infiltration) utilisée
ailleurs dans la littérature streambed. `ginette_velocity_for_row()` lit `vzm` tel quel ;
`darcy_velocity_from_head()` (repli analytique, convention opposée par construction) est negée
au point d'appel pour rester cohérente. Le régime d'identifiabilité (Pe<0.5/Pe>5) compare
`abs(peclet)`, le signe ne portant que la direction.

## Compatibilité Mac / Linux / cluster

- Chemins construits avec `pathlib`/`os.path`, aucun chemin absolu codé en dur dans le
  pipeline principal (0_ à 4_, `config_lomos.py`, `stallman_diffusivity.py`).
- `MAX_WORKERS` (`config_lomos.py`) plafonne le nombre de processus parallèles, calculé via
  `os.sched_getaffinity(0)` (respecte l'allocation cgroup sur un cluster Linux) avec repli sur
  `os.cpu_count()` (macOS, où `sched_getaffinity` n'existe pas).
- Le backend matplotlib bascule automatiquement en `Agg` (non interactif) uniquement sur Linux
  sans `$DISPLAY` (job cluster/SLURM en batch) - laissé intact sur macOS/Windows où l'absence de
  `$DISPLAY` ne signifie pas l'absence d'écran.
- `gfortran` doit être sur le `PATH` (`compile_ginette_src`, `src/src_python/Init_folders.py`)
  pour la compilation automatique du binaire `ginette` au premier lancement - installer via
  Homebrew sur Mac (`brew install gcc`) ou charger le module adéquat sur un cluster.