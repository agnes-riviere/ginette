"""
Script pour visualiser les résultats de simulation avec les fonctions de plotting du src.
Utilise les fonctions de Plot.py pour tracer les vitesses en 2D, les profils de température, flux, etc.
"""

import sys
import os
import pandas as pd
import numpy as np

# Ajouter le chemin vers src_python pour importer les fonctions de plotting
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/src_python'))

from Plot import (
    plot_water_profiles_interpol,
    plot_temperature_profiles_interpol,
    plot_heat_flux_profiles_interpolated,
    plot_fluxes_timeseries,
    plot_temperatures_profiles,
    plot_compare_temperatures_obs_sim
)

# ==================== Configuration ====================
# Définir la date de début de simulation
DATE_SIMUL_BEGIN = '2013-06-01'  # À adapter selon votre cas
FONTSIZE = 15

# ==================== Tracer les vitesses en 2D ====================
def plot_water_velocity_2d():
    """Trace les profils de vitesse d'eau en 2D (profondeur x temps)"""
    print("Traçage des profils de vitesse d'eau en 2D...")
    plot_water_profiles_interpol(fontsize=FONTSIZE, date_simul_bg=DATE_SIMUL_BEGIN)

# ==================== Tracer les profils de température en 2D ====================
def plot_temperature_2d():
    """Trace les profils de température en 2D (profondeur x temps)"""
    print("Traçage des profils de température en 2D...")
    plot_temperature_profiles_interpol(fontsize=FONTSIZE, date_simul_bg=DATE_SIMUL_BEGIN)

# ==================== Tracer les flux de chaleur en 2D ====================
def plot_heat_flux_2d():
    """Trace les profils de flux de chaleur en 2D (profondeur x temps)"""
    print("Traçage des profils de flux de chaleur en 2D...")
    plot_heat_flux_profiles_interpolated(fontsize=FONTSIZE, date_simul_bg=DATE_SIMUL_BEGIN)

# ==================== Tracer les séries temporelles des flux ====================
def plot_fluxes_ts():
    """Trace les séries temporelles des vitesses et flux à la profondeur max"""
    print("Traçage des séries temporelles des flux...")
    plot_fluxes_timeseries(fontsize=FONTSIZE, date_simul_bg=DATE_SIMUL_BEGIN)

# ==================== Tracer les profils de température 1D ====================
def plot_temperature_profiles_1d():
    """Trace l'évolution des profils de température à différents intervalles"""
    print("Traçage des profils de température 1D...")
    plot_temperatures_profiles(fontsize=FONTSIZE, interval=43200)  # Intervalle : 12 heures

# ==================== Main ====================
if __name__ == "__main__":
    print("=" * 60)
    print("VISUALISATION DES RÉSULTATS DE SIMULATION GINETTE")
    print("=" * 60)
    
    # Décommenter les fonctions qu'on souhaite exécuter :
    
    # Tracer les vitesses en 2D
    plot_water_velocity_2d()
    
    # Tracer les profils de température en 2D
    # plot_temperature_2d()
    
    # Tracer les flux de chaleur en 2D
    # plot_heat_flux_2d()
    
    # Tracer les séries temporelles
    # plot_fluxes_ts()
    
    # Tracer les profils 1D
    # plot_temperature_profiles_1d()
    
    print("=" * 60)
    print("Visualisations terminées !")
    print("=" * 60)