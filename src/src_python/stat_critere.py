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

def mse(predictions, ref,sigma):
    #Calculation of the differences between predictions and reference
    differences = ((predictions - ref)/sigma) 
    #Squares of differences                
    differences_squared = differences ** 2
    #sum of the differences squared                  
    sum_of_differences_squared = differences_squared.sum() 
    # # The root square of the mean
    # mse_val = np.sqrt(mean_of_differences_squared)   
    # Misfit value         
    return sum_of_differences_squared


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
