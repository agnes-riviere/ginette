# -*- coding: utf-8 -*-
"""
Grid Search Analysis for RIV2D Model
====================================
This script performs a grid search analysis to calibrate hydrological parameters
by comparing simulated and observed temperature data.

Author: Agnès Rivière
Date: 2023-10-01
Version: 1.0
Purpose: Parameter estimation using grid search methodology
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import importlib

# Configuration: Set up paths and station information
dir_ginette = "/home/ariviere/Programmes/ginette"
Station = "AmB"  # Example station, replace with actual station name
sys.path.append(dir_ginette)  # Add parent folder to Python path

# Import required modules from the ginette package
from src.src_gmsh import mesh_generator
from src.src_python import Init_folders
from src.src_python import Direct_model
from src.src_python import Read_obs
from src.src_python import Plot
from src.src_python import stat_critere
from src.src_python import Grid_search

# Import all functions from modules (use with caution in production)
from src.src_python.Init_folders import *
from src.src_python.Direct_model import *
from src.src_python.Read_obs import *
from src.src_python.Plot import *
from src.src_python.stat_critere import *
from src.src_python.Grid_search import *

# Reload modules for development (ensures latest changes are loaded)
importlib.reload(Init_folders)
importlib.reload(Direct_model)
importlib.reload(Grid_search)

print('Working directory:', dir_ginette)

# Step 1: Read configuration from input file
# ==========================================
main_dir = os.path.join(dir_ginette, "application/RIV2D", Station)
os.chdir(main_dir)
print("Current directory:", main_dir)

# Read and parse the input configuration file
with open("input_inversion.txt", "r") as f:
    lines = f.readlines()

def parse_line(line):
    """Helper function to parse configuration lines of the form 'key: value'"""
    key, value = line.strip().split(":", 1)
    return key.strip(), value.strip()

# Parse all configuration parameters
parsed = dict(parse_line(line) for line in lines)

# Extract and convert configuration parameters
date_simul_bg = parsed["date_simul_bg"]
Station = parsed["Station"]
param_struct = [x.strip(" '") for x in parsed["param_struct"].split(",")]
sensors = [x.strip(" '") for x in parsed["sensors"].split(",")]
sigma = float(parsed["sigma"])
zones_to_invert = [int(z) for z in parsed["zones_to_invert"].strip("[]").split(",")]
parameters_to_invert = [x.strip(" '") for x in parsed["parameters_to_invert"].split(",")]
simul_todo_range = parsed["simul_todo"].split(",")
simul_todo = range(int(simul_todo_range[0]), int(simul_todo_range[1]))
main_dir = parsed["path_simul"]
dir_ginette = parsed["dir_ginette"]

os.chdir(main_dir)
main_dir = os.path.join(dir_ginette, "application/RIV2D", Station)

# Convert date string to datetime object for time series analysis
date_simul_bg = pd.to_datetime(date_simul_bg, format='%Y-%m-%d %H:%M:%S')

# Display parsed configuration for verification
print(f"Simulation start date: {date_simul_bg}")
print(f"Station: {Station}")
print(f"Parameter structure: {param_struct}")
print(f"Sensors: {sensors}")
print(f"Sigma: {sigma}")
print(f"Zones to invert: {zones_to_invert}")
print(f"Parameters to invert: {parameters_to_invert}")
print(f"Simulations to process: {simul_todo}")
print(f"Simulation path: {main_dir}")
print(f"Ginette directory: {dir_ginette}")

# Step 2: Load input data
# =======================

# Load parameter table containing all parameter combinations to test
param_table = pd.read_csv("param_table.txt", sep=",", header=0)

# Load observed temperature data from sensors
obs_temp = pd.read_csv("Obs_temp_PT100_t.dat", sep=",", header=0)

# Clean and prepare observed temperature data
for i in range(len(sensors)):
    # Convert sensor data to numeric, replacing invalid values with NaN
    obs_temp[sensors[i]] = pd.to_numeric(obs_temp[sensors[i]], errors='coerce')

# Ensure time column is numeric (in seconds)
obs_temp['Time'] = pd.to_numeric(obs_temp['Time'], errors='coerce')

# Convert dates to datetime objects for proper time series handling
obs_temp['dates'] = pd.to_datetime(obs_temp['dates'], format='%Y-%m-%d %H:%M:%S')

# Remove initial equilibrium period (first two days) from observations
obs_temp = remove_first_two_days_obs(obs_temp, date_simul_bg)

# Step 3: Initialize result storage
# =================================

# Add simulation index to parameter table
param_table['index_sim'] = param_table.index + 1
header_save_param = param_table.columns

# Initialize tables to store results
param_table_simul = pd.DataFrame(columns=header_save_param)
mse_table = pd.DataFrame(columns=['index_sim'] + sensors + ['Total_mse'])

# Step 4: Check for missing simulation files
# ===========================================
missing_simulations = []
for id_sim in simul_todo:
    output_file = f"OUTPUT_{Station}/{id_sim}_S_temp_PT100_t.dat"
    if not os.path.exists(output_file):
        missing_simulations.append(id_sim)

if missing_simulations:
    print(f"Missing simulations: {missing_simulations}")
    sys.exit("Script stopped due to missing simulations.")

# Step 5: Main analysis loop - Grid search execution
# ==================================================
print("Starting grid search analysis...")

for id_sim in simul_todo:
    # Analyze each simulation: compute MSE between observed and simulated data
    mse, param_table_id_sim = analysis_gridsearch_2D(
        Station, main_dir, date_simul_bg, obs_temp, sensors,
        param_table, id_sim, param_struct, zones_to_invert,
        parameters_to_invert, param_table_simul, sigma
    )
    
    # Store results: MSE values and corresponding parameters
    mse['index_sim'] = id_sim
    mse_table = pd.concat([mse_table, mse], ignore_index=True)
    param_table_simul = pd.concat([param_table_simul, param_table_id_sim], ignore_index=True)

# Step 6: Combine MSE results with parameter values
# =================================================
for param in header_save_param:
    mse_table[param] = param_table_simul[param].values

# Step 7: Statistical analysis and normalization
# ==============================================

# Calculate normalization factors
nb_it_time = obs_temp['Time'].count()  # Number of time steps
nb_sensors = len(sensors)              # Number of sensors

print(f"Number of time iterations: {nb_it_time}")
print(f"Number of sensors: {nb_sensors}")

# Normalize total MSE by number of observations
normalize_total_mse(mse_table, col="Total_mse", 
                   normalized_col="Total_mse_normalized",
                   nb_it_time=nb_it_time*nb_sensors)

# Normalize individual sensor MSE values
for sensor in sensors:
    normalize_total_mse(mse_table, col=f"{sensor}", 
                       normalized_col=f"{sensor}_normalize_mse",
                       nb_it_time=nb_it_time)

# Step 8: Save results
# ===================

# Save complete results table
mse_table.to_csv("SENSI_" + Station + "/" + "S_mse_simul.dat", index=False)
print("Results saved to:", "SENSI_" + Station + "/" + "S_mse_simul.dat")
print(mse_table.head())

# Step 9: Identify best model (lowest total MSE)
# ==============================================
best_row = mse_table.loc[mse_table['Total_mse'].idxmin()]

# Extract parameters of interest for the best model
params_of_interest = best_row[["k4", "n4", "l4", "k5", "n5", "l5", 'Total_mse']]
print("Parameter values for the best model (lowest Total_mse):")
print(params_of_interest)

# Save best parameters
params_of_interest.to_csv("SENSI_" + Station + "/" + "S_best.dat", index=False)

print("Analysis completed successfully!")
