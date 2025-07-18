# -*- coding: utf-8 -*-
"""
This module provides statistical criteria and metrics for comparing simulated and observed data,
particularly for temperature sensors. It includes functions to compute RMSE, MSE, and a summary
of metrics (RMSE, MAE, PBias, KGE) for multiple sensors.
Functions:
----------
rmse(predictions, ref, sigma):
    Computes the Root Mean Square Error (RMSE) between predictions and reference values,
    normalized by sigma.
mse(predictions, ref, sigma):
    Computes the Mean Squared Error (MSE) between predictions and reference values,
    normalized by sigma.
calculate_metrics(data):
    Calculates several statistical metrics (RMSE, MAE, PBias, KGE) for three temperature sensors
    ('Temp1', 'Temp2', 'Temp3') based on simulated and observed data columns in the input DataFrame.
Dependencies:
-------------
- numpy
- pandas
- collections.namedtuple
- dataclasses.dataclass
- random.uniform, random.gauss
- typing
"""
import numpy as np
import pandas as pd
from collections import namedtuple
from dataclasses import dataclass
from random import uniform, gauss
from numpy import inf, nansum, log, size, var, mean, isclose, sqrt, zeros, all
from typing import Callable, Tuple, Sequence

def rmse(predictions, ref,sigma):
    differences = ((predictions - ref)/sigma)                      #the DIFFERENCEs.
    differences_squared = differences ** 2                    #the SQUAREs of ^
    mean_of_differences_squared = differences_squared.mean()  #the MEAN of ^
    rmse_val = np.sqrt(mean_of_differences_squared)           #ROOT of ^

    return rmse_val
#_____________________________________________________
import numpy as np
import pandas as pd

def mse(predictions, ref, sigma):
    """
    Calcule la somme des carrés normalisés des erreurs (misfit),
    souvent utilisée pour la vraisemblance sous hypothèse gaussienne.

    Arguments :
    - predictions : np.ndarray ou pd.Series des valeurs prédites
    - ref : np.ndarray ou pd.Series des valeurs de référence (observées)
    - sigma : float ou np.ndarray (écart type du bruit)

    Retourne :
    - float : somme des erreurs quadratiques normalisées
    """
    # Calcul des différences normalisées
    differences = (predictions - ref) / sigma
    # Carré des différences
    differences_squared = differences ** 2
    # Somme des carrés
    sum_of_differences_squared = differences_squared.sum()
    return sum_of_differences_squared
#_____________________________________________________
def normalize_total_mse(mse_table, col="Total_mse", normalized_col="Total_mse_normalized", nb_it_time=2688):
    """
    Ajoute une colonne à mse_table contenant la colonne total_col normalisée par divisor.

    Arguments :
    - mse_table : pd.DataFrame
    - total_col : str, nom de la colonne à normaliser
    - normalized_col : str, nom de la nouvelle colonne normalisée
    - divisor : float, valeur de normalisation

    Retourne :
    - pd.DataFrame avec la nouvelle colonne ajoutée
    """
    mse_table[normalized_col] = mse_table[col] / nb_it_time
    return mse_table
#_____________________________________________________

def likelihood(mse_table, col_name, add_column=True):
    """
    Calcule la vraisemblance sous l’hypothèse que chaque MSE suit une distribution normale.

    Arguments :
    - mse_table : pd.DataFrame contenant les MSE
    - col_name : str, nom de la colonne contenant les MSE
    - add_column : bool, si True ajoute une colonne avec les vraisemblances

    Retourne :
    - np.ndarray : tableau des vraisemblances
    """
    if col_name not in mse_table.columns:
        raise KeyError(f"Column '{col_name}' not found in the DataFrame.")

    name_col_likelihood = f"{col_name}_likelihood"
    mse_values = mse_table[col_name].astype(np.float64)

    # Vraisemblance gaussienne sous forme non normalisée (log-vraisemblance simplifiée)
    likelihoods = np.exp(-0.5 * mse_values)

    if add_column:
        mse_table[name_col_likelihood] = likelihoods

    return likelihoods
#_____________________________________________________
def log_likelihood(mse_table, col_name, add_column=True):
    """
    Calcule la log-vraisemblance à partir de MSE.
    logL = -0.5 * MSE (forme simplifiée sans constante)
    """
    if col_name not in mse_table.columns:
        raise KeyError(f"Column '{col_name}' not found in the DataFrame.")
    
    log_likes = -0.5 * mse_table[col_name].astype(np.float64)
    
    if add_column:
        mse_table[f"{col_name}_log_likelihood"] = log_likes
    
    return log_likes

#______________________________________________________
def calculate_metrics(data):
    metrics = {}

    for i in range(1, 4):
        sensor = f'Temp{i}'
        diff = data[f'{sensor}_sim'] - data[f'{sensor}_obs']
        
        rmse = np.sqrt(np.mean(diff**2))  # RMSE
        mae = np.mean(np.abs(diff))  # MAE
        bias = 100 * np.sum(diff) / np.sum(np.abs(data[f'{sensor}_obs']))  # PBias

        std_obs = data[f'{sensor}_obs'].std()
        std_sim = data[f'{sensor}_sim'].std()
        mean_obs = data[f'{sensor}_obs'].mean()
        mean_sim = data[f'{sensor}_sim'].mean()

        kge = 1 - np.sqrt((np.square(std_sim/std_obs - 1)) +
                            (np.square(mean_sim/mean_obs - 1)) +
                            (np.square(rmse/std_obs - 1)))

        metrics[sensor] = {'RMSE': rmse, 'MAE': mae, 'PBias': bias, 'KGE': kge}

    return metrics


def calcul_indicateurs(df, t_event, var_forcage, var_reponse, window='2D'):
    """
    Calcule des indicateurs de réactivité autour d'un événement
    - df : dataframe avec index datetime
    - t_event : timestamp de l'événement (pluie, canicule)
    - var_forcage : ex. 'pluie' ou 'temp_air'
    - var_reponse : ex. 'hauteur' ou 'temp_riviere'
    - window : durée autour de l'événement à analyser
    """
    delta = pd.Timedelta(window)
    subset = df.loc[t_event - delta : t_event + delta]

    if subset.empty or subset[var_reponse].isna().all():
        return None

    # 1. Valeur initiale et pic de réponse
    val_ini = subset.loc[:t_event, var_reponse].iloc[-1]
    val_max = subset[var_reponse].max()
    val_min = subset[var_reponse].min()
    t_max = subset[var_reponse].idxmax()
    t_min = subset[var_reponse].idxmin()

    # 2. Délai de réponse
    lag_max = (t_max - t_event).total_seconds() / 3600
    lag_min = (t_min - t_event).total_seconds() / 3600

    # 3. Amplitude
    amp_max = val_max - val_ini
    amp_min = val_min - val_ini

    # 4. Temps de retour à ±10%
    seuil_sup = val_ini + 0.1 * amp_max if amp_max > 0 else val_ini + 0.1 * amp_min
    seuil_inf = val_ini - 0.1 * abs(amp_min)

    retour = subset.loc[t_event:]
    t_retour = retour[(retour[var_reponse] < seuil_sup) & (retour[var_reponse] > seuil_inf)]

    if not t_retour.empty:
        retour_time = (t_retour.index[0] - t_event).total_seconds() / 3600
    else:
        retour_time = np.nan

    # 5. Taux de variation
    duration_h = (subset.index[-1] - subset.index[0]).total_seconds() / 3600
    slope = (subset[var_reponse].iloc[-1] - val_ini) / duration_h

    # 6. Aire sous la courbe (effet cumulatif)
    area = np.trapz(y=subset[var_reponse] - val_ini, x=(subset.index - subset.index[0]).total_seconds() / 3600)

    return {
        'event_time': t_event,
        'lag_max_hr': lag_max,
        'lag_min_hr': lag_min,
        'amp_max': amp_max,
        'amp_min': amp_min,
        'retour_hr': retour_time,
        'slope_hr': slope,
        'area_integrated': area
    }
