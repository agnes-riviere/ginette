#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures de résultats du grid search pour un cas réel (lomos230) : marginales 1D,
cartes 2D du misfit par paire de paramètres, et comparaison temporelle
température observée / simulée pour la meilleure combinaison de paramètres.

Le nombre de paramètres calibrés (Name_parameters) est variable : le nombre de
marginales et de cartes 2D (une par paire de paramètres) s'adapte automatiquement,
comme la construction de la grille dans 1_define_grid_search.py.

Contrairement au cas synthétique (SYNTHETIC_CASES), il n'y a pas de "vraie
valeur" connue pour les paramètres dans un cas réel : seule la meilleure
combinaison trouvée (misfit minimal) est mise en évidence.
"""

# %% IMPORTS
import sys
from itertools import combinations
from pathlib import Path

project_root = Path(__file__).resolve().parents[2]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

cmap = matplotlib.colormaps['nipy_spectral'].resampled(256).reversed()
plt.rcParams["font.size"] = 12

BASE_APP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE_APP_DIR))
# POINT_NAME partagé via config_lomos.py avec les autres scripts (0_/2_/3_) -
# avant ce fix, ce script gardait son propre POINT_NAME="lomos231" en dur,
# désynchronisé de celui utilisé pour générer les résultats ("lomos230").
from config_lomos import POINT_NAME
RESULTS_DIR = BASE_APP_DIR / "results" / POINT_NAME

# %% PARAMÈTRES UTILISATEUR

# Doit correspondre à Name_parameters de 1_define_grid_search.py (même liste,
# même ordre pas nécessaire). Le nombre de paramètres peut varier librement :
# le nombre de marginales et de cartes 2D s'adapte automatiquement ci-dessous.
Name_parameters = ["log_k", "lam", "n"]

# Deux critères de misfit sont calculés par 3_misfit.py (voir son en-tête) :
# misfit_tot (quadratique/RMSE, sensible aux petits décalages de phase) et
# mae_tot (erreur absolue moyenne, plus robuste à ces décalages - cf. Cognac &
# Ronayne 2023). Ils NE désignent pas toujours la même combinaison comme
# meilleure : on produit donc les figures pour LES DEUX, avec un suffixe dans
# le nom de fichier, plutôt que de choisir arbitrairement l'un des deux.
ASSESS_VARS = ["misfit_tot", "mae_tot"]

# %% CHARGEMENT DES RÉSULTATS
results = pd.read_csv(RESULTS_DIR / "results.txt", sep=" ", index_col=[0])
obs_data = pd.read_csv(RESULTS_DIR / "observed_data.txt", sep=" ", index_col=[0])
pairs = list(combinations(Name_parameters, 2))

# Cache de toutes les simulations dont le misfit a pu être calculé (chargées
# une seule fois, réutilisées pour les deux critères de misfit - évite de
# relire tous les fichiers sim_temp_*.txt deux fois).
valid_ids = results.loc[results[ASSESS_VARS[0]].notna(), "ID"].astype(int).tolist()
sim_cache = {sid: pd.read_csv(RESULTS_DIR / f"sim_temp_{sid}.txt", sep=" ", index_col=[0]) for sid in valid_ids}
n_missing = len(results) - len(valid_ids)
if n_missing:
    print(f"Attention : {n_missing} simulations sur {len(results)} n'ont pas de misfit calculé "
          "(fichier manquant/vide) - absentes des figures ci-dessous.")

bests = {}
for assess_var in ASSESS_VARS:
    best = results.loc[results[assess_var].idxmin(), "ID"]
    bests[assess_var] = best
    print(f"Meilleure combinaison selon {assess_var} "
          f"(ID={int(best)}, {assess_var}={results.loc[results['ID'] == best, assess_var].values[0]:.1f}) :")
    for param in Name_parameters:
        print(f"  {param} = {results.loc[results['ID'] == best, param].values[0]}")

if bests[ASSESS_VARS[0]] != bests[ASSESS_VARS[1]]:
    print(f"\nATTENTION : {ASSESS_VARS[0]} et {ASSESS_VARS[1]} désignent des combinaisons "
          "DIFFÉRENTES comme meilleures - voir les figures des deux critères ci-dessous.")

# %% FIGURES INDIVIDUELLES (une figure par marginale, par carte 2D, et par
# profondeur pour la série temporelle - plus simple à consulter/exporter
# séparément qu'un seul grand mosaïque), répétées pour CHAQUE critère de
# misfit. Chaque figure est sauvegardée sous son propre nom dans results/ ET
# affichée avant de passer à la suivante.
for assess_var in ASSESS_VARS:
    best = bests[assess_var]
    sim_data = pd.read_csv(RESULTS_DIR / f"sim_temp_{int(best)}.txt", sep=" ", index_col=[0])

    # Marginales 1D : somme du critère pour chaque valeur testée d'un paramètre
    # (une valeur avec un faible cumul marginal est globalement mieux supportée
    # par les données, tous les autres paramètres confondus)
    marginals = {param: results.groupby(param)[assess_var].sum() for param in Name_parameters}

    # --- Marginales 1D (une figure par paramètre) ---
    for param in Name_parameters:
        fig, ax = plt.subplots(figsize=(5, 4), dpi=110)
        m = marginals[param]
        ax.plot(m.index, m.values, 'o-', lw=2, color="black")
        ax.axvline(results.loc[results['ID'] == best, param].values[0], lw=2, ls="--",
                   color="gold", label="Meilleure valeur")
        ax.set_xlabel(param)
        ax.set_ylabel(f"{assess_var} (marginale)")
        ax.set_title(f"Marginale : {param} ({assess_var})")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(RESULTS_DIR / f"marginal_{param}_{assess_var}.png", dpi=150)
        plt.show()

    # --- Cartes 2D du misfit (une figure par paire de paramètres) ---
    for p1, p2 in pairs:
        fig, ax = plt.subplots(figsize=(5.5, 4.5), dpi=110)
        sca = ax.scatter(results[p1], results[p2], c=results[assess_var], cmap=cmap)
        fig.colorbar(sca, ax=ax, label=assess_var)
        ax.scatter(results.loc[results['ID'] == best, p1], results.loc[results['ID'] == best, p2],
                   marker="*", edgecolors="black", facecolors="gold", s=150, lw=1, label="Meilleure combinaison")
        ax.set_xlabel(p1)
        ax.set_ylabel(p2)
        ax.set_title(f"{assess_var} : {p1} vs {p2}")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(RESULTS_DIR / f"misfit2d_{p1}_{p2}_{assess_var}.png", dpi=150)
        plt.show()

    # --- Série temporelle : observations vs meilleure simulation (une figure par profondeur) ---
    # obs_data couvre toute la période disponible (observed_data.txt), qui s'étend souvent
    # bien au-delà de la durée simulée (nb_day) : on restreint donc l'affichage à la
    # période RÉELLEMENT simulée, sinon la courbe observée qui continue seule au-delà de
    # la fin de la simulation donne l'illusion trompeuse d'une dérive/divergence du modèle.
    t_max_sim_days = sim_data.Time.max() / 86400
    obs_plot = obs_data[obs_data.Time / 86400 <= t_max_sim_days]

    temp_cols = [c for c in ["Temp1", "Temp2", "Temp3"] if c in obs_plot.columns and c in sim_data.columns]
    for col in temp_cols:
        fig, ax = plt.subplots(figsize=(10, 4), dpi=110)
        ax.plot(obs_plot.Time / 86400, obs_plot[col], ls="-", lw=2, color="tab:blue", label=f"Observé {col}")
        ax.plot(sim_data.Time / 86400, sim_data[col], ls="--", lw=1.5, color="tab:red",
                label=f"Simulé {col} (best ID={int(best)}, {assess_var})")
        ax.set_xlabel("Temps depuis le début de la simulation [jours]")
        ax.set_ylabel("Température [°C]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
        ax.set_title(f"Meilleure simulation vs observations : {col} ({POINT_NAME}, {assess_var})")
        fig.tight_layout()
        fig.savefig(RESULTS_DIR / f"timeseries_{col}_best_{assess_var}.png", dpi=150)
        plt.show()

    # --- Toutes les simulations superposées, colorées par le misfit (une
    # figure, un sous-graphique par profondeur) : donne une vue d'ensemble de
    # comment l'ensemble du balayage se comporte par rapport aux observations,
    # pas seulement le meilleur point. Tracées du pire au meilleur misfit pour
    # que les meilleurs ajustements restent visibles au premier plan.
    norm = matplotlib.colors.Normalize(vmin=results[assess_var].min(), vmax=results[assess_var].max())
    order = (results.loc[results['ID'].isin(sim_cache), ['ID', assess_var]]
             .sort_values(assess_var, ascending=False)['ID'].astype(int).tolist())

    fig, axs = plt.subplots(len(temp_cols), 1, figsize=(11, 3.2 * len(temp_cols)), sharex=True, dpi=110)
    if len(temp_cols) == 1:
        axs = [axs]
    for sid in order:
        sim_i = sim_cache[sid]
        color = cmap(norm(results.loc[results['ID'] == sid, assess_var].values[0]))
        for ax, col in zip(axs, temp_cols):
            ax.plot(sim_i.Time / 86400, sim_i[col], color=color, lw=0.6, alpha=0.7)
    for ax, col in zip(axs, temp_cols):
        ax.plot(obs_plot.Time / 86400, obs_plot[col], color="black", lw=2.2, label="Observé", zorder=10)
        ax.set_ylabel("Température [°C]")
        ax.set_title(col)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9, loc="upper left")
    axs[-1].set_xlabel("Temps depuis le début de la simulation [jours]")
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    fig.colorbar(sm, ax=axs, label=assess_var, shrink=0.8)
    fig.suptitle(f"Toutes les simulations ({len(order)}) colorées par {assess_var}")
    fig.savefig(RESULTS_DIR / f"spaghetti_{assess_var}.png", dpi=150)
    plt.show()
