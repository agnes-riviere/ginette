CALAGE DE GINETTE SUR DONNEES REELLES LOMOS-mini (lomos230)
=============================================================

Ce dossier calibre Ginette sur des mesures de terrain (dispositif LOMOS-mini, Cucchi et al.
2018 : 4 sondes de temperature a 10/20/30/40cm de profondeur dans le lit d'une riviere, un
capteur de temperature en eau libre TempMolo, et un capteur de difference de charge deltaP)
pour estimer les echanges nappe-riviere (permeabilite, porosite, conductivite thermique du
streambed).


PIPELINE
--------

1. 0_boundary_conditions_real_case.py - lit le CSV d'observation, construit les conditions
   initiales et aux limites, ecrit les fichiers d'entree Ginette (GINETTE_SENSI/).
2. 1_define_grid_search.py - genere results/grid_search.csv, la grille de parametres a
   tester (log_k, lam, n, c).
3. 2_run_real_case.py - lance une simulation Ginette par ligne de la grille (parallelise,
   MAX_WORKERS=6 pour ne pas saturer un poste de travail local).
4. 3_misfit.py - compare chaque simulation aux observations, calcule les metriques
   d'ajustement (misfit L2, MAE, KGE, pbias) et un diagnostic d'identifiabilite (nombre de
   Peclet, largeur du plateau de misfit).
5. 4_plot_results.py - figures de synthese.

Tous les parametres utilisateur (point d'observation, config BC A/C, periode de simulation,
MAX_WORKERS, SAVE_VELOCITY_PROFILES, bornes et pas du grid search) sont centralises dans
config_lomos.py - un seul fichier a editer, partage par les 5 scripts du pipeline et
stallman_diffusivity.py.

La periode de simulation (DATE_SIMUL_BG/NB_DAY) peut etre imposee a la main
(FORCE_DATE_SIMUL_BG=True, valeurs ci-dessous) ou reprise automatiquement de la
recommandation de stallman_diffusivity.py (FORCE_DATE_SIMUL_BG=False : relit
results/stallman_recommended_period.txt ; si ce fichier n'existe pas encore, config_lomos.py
lance lui-meme stallman_diffusivity.py pour le generer).

stallman_diffusivity.py lit POINT_NAME depuis config_lomos.py comme les autres scripts et
marche donc sur n'importe quel point au format LOMOS-mini (lomos230, lomos231, ...) ; seules
les profondeurs de capteurs (10/20/30/40cm) restent codees en dur, propres au dispositif
LOMOS-mini.

config_lomos.py controle aussi l'enchainement a partir de 2_run_real_case.py :
  - RUN_GRID_SEARCH : relance les simulations Ginette (False = reutilise les sim_temp_*.txt
    deja presents, utile pour refaire juste le misfit/les plots).
  - RUN_MISFIT : enchaine sur 3_misfit.py apres le grid search.
  - RUN_PLOTS : enchaine sur 4_plot_results.py apres le misfit.

Chaque script ne declenche que l'etape suivante (2->3->4), donc lancer 3_misfit.py seul avec
RUN_PLOTS=True regenere aussi les figures, sans repasser par 2_run_real_case.py.


LE PROBLEME DE LA CONDITION LIMITE HAUTE : TempMolo n'est pas la temperature du sediment
------------------------------------------------------------------------------------------

TempMolo est mesure 1-2cm au-dessus du lit, en eau libre - ce n'est pas la meme grandeur
physique que la temperature du milieu poreux a z=0 (melange convectif de l'eau libre vs
conduction/advection dans le sediment). L'utiliser tel quel comme condition de Dirichlet en
surface revient a imposer au modele un forcage dont il ne peut pas reproduire correctement la
physique, quels que soient les parametres testes.

Le choix des capteurs de CL vit dans config_lomos.py (TOP_SENSOR/BOTTOM_SENSOR, librement
choisis parmi SENSOR_DEPTHS) : domaine, points de comparaison, rescaling de deltaP en
decoulent tous automatiquement (voir la section Pipeline plus haut). Deux presets deja
etudies sur lomos230 :

  - Config A (TOP_SENSOR='TempMolo', BOTTOM_SENSOR='Temp4') : TempMolo brut en
    CL haute, Temp4 en CL basse, domaine 0 a -0.4m, calage sur Temp1/Temp2/Temp3.
  - Config C (TOP_SENSOR='Temp1', BOTTOM_SENSOR='Temp4', par defaut, recommandee) : Temp1 ET Temp4, mesures reelles
    sans aucune correction, en CL haute/basse, domaine 0 a -0.3m, calage sur Temp2/Temp3
    seulement. deltaP est ramene au sous-domaine (x0.75) en supposant un gradient
    hydraulique homogene sur la colonne.

Le nombre de points d'observation Ginette reste fige a 4, equidistants de DZ_OBS sous
TOP_SENSOR (format E_parametre_bck.dat, partage avec d'autres apps) : DZ_OBS est derive de
SENSOR_DEPTHS, qui doit donc rester uniformement espace (ex: 10cm entre capteurs pour
LOMOS-mini) - seule la profondeur/l'espacement peuvent varier d'une station a l'autre, pas le
nombre de capteurs.

Une troisieme option (TempMolo corrige par extrapolation Stallman depuis la paire fiable
Temp1-Temp2) a ete testee puis abandonnee : la correction ne touche que l'harmonique 24h,
laissant la tendance sur plusieurs jours (dominante sur 42 jours de calage) inchangee.
Etendre la correction a tout le spectre par FFT est mal pose - le gain de Stallman croit avec
la racine de la frequence et amplifie le bruit haute frequence de facon exponentielle,
produisant des temperatures non physiques meme avec un filtre passe-bas.


CONFIG A ET CONFIG C NE SONT PAS JUSTE "DIFFERENTES" : ELLES N'IDENTIFIENT PAS LA MEME CHOSE
------------------------------------------------------------------------------------------------

Le nombre de Peclet thermique (Pe = flux advectif / flux conductif, Kurylyk et al. 2017
eq.12) dit a priori quels parametres sont identifiables a partir de la seule temperature
(Cucchi, Flipo, Riviere, Rubin 2021, qui utilise Ginette sur ce meme type de site) :
  - Pe < 0.5 (conductif) : seule la diffusivite thermique effective est identifiable, la
    permeabilite seulement bornee.
  - 0.5 < Pe < 5 (transition) : permeabilite ET diffusivite identifiables.
  - Pe > 5 (advectif) : permeabilite identifiable, diffusivite non identifiable.

Sur lomos230 (42 jours, 2026-04-30 a 2026-06-11) :

  Config A : log_k best-fit = -16, |Pe| ~ 2e-4, plateau top-5% du misfit sur log_k =
             [-18, -14] (4 ordres de grandeur)
             -> log_k NON identifiable - "-16" est un point pris au hasard sur un plateau
                plat.
  Config C : log_k best-fit = -12, |Pe| ~ 0.63, plateau top-5% du misfit sur log_k = {-12}
             (une seule valeur de grille)
             -> log_k bien identifie.

Config C fit aussi mieux sur les capteurs communs, en comparaison honnete (le kge_3 de C
compare la maille de bord a sa propre condition limite, donc trivial, a exclure) : KGE
Temp2/Temp3 = 0.94/0.96 pour C contre 0.84/0.93 pour les memes profondeurs en A. Et le
log_k=-12 de Config C est coherent avec l'estimation Stallman independante (analyse
spectrale pure, sans Ginette) : -13.4 a -12.1 (stallman_diffusivity.py).

  NB : chiffres recalcules apres les corrections du 2026-07-24 (voir section suivante) ; a
  reconfirmer si tu relances un grid search complet, notamment log_k best-fit qui depend
  indirectement de Cv via le spin-up automatique.

Config C est recommandee par defaut sur ce site. Config A reste disponible mais son log_k
doit etre lu comme une borne, pas une estimation. lam et n restent mal contraints
individuellement meme en Config C (equifinalite) - seule leur combinaison (diffusivite
effective) semble contrainte par les donnees ; eviter de sur-interpreter leurs valeurs
ponctuelles.


GRILLE DE PARAMETRES
---------------------

Volontairement large et non-informative meme la ou on a une idee a priori des parametres
(log_k sur 7 ordres de grandeur : -18 a -11), pour pouvoir caracteriser l'incertitude/le
plateau de misfit plutot que de la restreindre a priori autour d'une estimation.

c (capacite calorifique SPECIFIQUE du solide, J/kg/K, fixe a 2000 - PAS volumique malgre ce
que suggere son nom court) est combinee a rhosi (densite du solide, E_p_therm.dat, fixe a
1180 kg/m3) via CASE('ZHZ') (src/ginette_V2.f90). 1180 kg/m3 (baisse depuis 2650/quartz pur
le 2026-07-24) reflete un sediment riche en matiere organique - coherent avec c=2000, qui est
eleve pour un solide purement mineral. Ensemble, melanges a l'eau selon la porosite
(Cv = n*Cw_vol + (1-n)*rhosi*c), ils donnent une capacite volumique du milieu sature de 2.45
a 3.46 MJ/m3/K sur n=[0.05,0.6], coherente avec la litterature streambed (~2-3.5 MJ/m3/K,
Stallman/Lapham/Constantz).


CONVENTION DE SIGNE DU FLUX (darcy_flux_m_s, peclet dans 3_misfit.py)
------------------------------------------------------------------------

Positif = vers le HAUT (exfiltration), negatif = vers le BAS (infiltration) - convention de
Ginette elle-meme (z positif vers le haut, src/ginette_V2.f90 ~l.3712 :
v = -K/mu*(dP/dz+rho*g)), pas celle de Stallman (z=profondeur, positif=infiltration) utilisee
ailleurs dans la litterature streambed. ginette_velocity_for_row() lit vzm tel quel ;
darcy_velocity_from_head() (repli analytique, convention opposee par construction) est negee
au point d'appel pour rester coherente. Le regime d'identifiabilite (Pe<0.5/Pe>5) compare
abs(peclet), le signe ne portant que la direction.


COMPATIBILITE MAC / LINUX / CLUSTER
--------------------------------------

  - Chemins construits avec pathlib/os.path, aucun chemin absolu code en dur dans le
    pipeline principal (0_ a 4_, config_lomos.py, stallman_diffusivity.py).
  - MAX_WORKERS (config_lomos.py) plafonne le nombre de processus paralleles, calcule via
    os.sched_getaffinity(0) (respecte l'allocation cgroup sur un cluster Linux) avec repli
    sur os.cpu_count() (macOS, ou sched_getaffinity n'existe pas).
  - Le backend matplotlib bascule automatiquement en Agg (non interactif) uniquement sur
    Linux sans $DISPLAY (job cluster/SLURM en batch) - laisse intact sur macOS/Windows ou
    l'absence de $DISPLAY ne signifie pas l'absence d'ecran.
  - gfortran doit etre sur le PATH (compile_ginette_src, src/src_python/Init_folders.py)
    pour la compilation automatique du binaire ginette au premier lancement - installer via
    Homebrew sur Mac (brew install gcc) ou charger le module adequat sur un cluster.
