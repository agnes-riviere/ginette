#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 09:51:06 2025

@author:  Agnès Rivière, Samuel Larance
"""


# IMPORT:
import sys
from pathlib import Path
# Add project root to path
project_root = Path(__file__).resolve().parents[2]
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))
import importlib
import os
import numpy as np
import pandas as pd
from time import time
import shutil
import multiprocessing as mp


# fonction plot resultats 3 temperature in the same plot of sim_temp and obs_temp colones 'Time' and 3 temperatures 'Temp1','Temp2','Temp3'
def plot_time_series(sim_temp, obs_temp, ID_simul):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    plt.plot(sim_temp.index, sim_temp['Temp1'], label='Simulated Temp1', color='blue')
    plt.plot(sim_temp.index, sim_temp['Temp2'], label='Simulated Temp2', color='orange')
    plt.plot(sim_temp.index, sim_temp['Temp3'], label='Simulated Temp3', color='green')
    plt.scatter(obs_temp.index, obs_temp['Temp1'], label='Observed Temp1', color='blue', marker='x')    
    plt.scatter(obs_temp.index, obs_temp['Temp2'], label='Observed Temp2', color='orange', marker='x')
    plt.scatter(obs_temp.index, obs_temp['Temp3'], label='Observed Temp3', color='green', marker='x')
    plt.xlabel('Time')
    plt.ylabel('Temperature (°C)')
    plt.title(f'Temperature Time Series for Simulation ID {ID_simul}')
    plt.legend()
    plt.grid()
    plt.show()

#observed_data.txt in results directory
obs_temp=pd.read_csv('results/observed_data.txt',sep=' ',index_col=0) 
sim_temp=pd.read_csv('results/sim_temp_0.txt',sep=' ',index_col=0) 

print('Observed Data:')
print(obs_temp.head())
plot_time_series(sim_temp, obs_temp, ID_simul=0)

