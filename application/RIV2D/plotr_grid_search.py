# -*- coding: utf-8 -*-
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import importlib

# Configuration: Set the main directory and station
dir_ginette = "/home/ariviere/Programmes/ginette"
Station = "AmB"  # Example station, replace with actual station name
sys.path.append(dir_ginette)  # Add parent folder to Python path

# Import RIV2D modules
from src.src_gmsh import mesh_generator
from src.src_python import Init_folders, Direct_model, Read_obs, Plot, stat_critere, Grid_search

# Import all functions/classes from the relevant modules
from src.src_gmsh.mesh_generator import *
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

print('Repository:', dir_ginette)

# Set working directory and read input parameters
main_dir = os.path.join(dir_ginette, "application/RIV2D", Station)
os.chdir(main_dir)
print('Working directory:', main_dir)

# Read and parse input configuration file
with open("input_inversion.txt", "r") as f:
    lines = f.readlines()

# Helper function to parse configuration lines
def parse_line(line):
    """Parse a line of the form 'key: value' and return key, value tuple"""
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

# Update working directory based on parsed configuration
os.chdir(main_dir)
main_dir = os.path.join(dir_ginette, "application/RIV2D", Station)

# Convert date to datetime format
date_simul_bg = pd.to_datetime(date_simul_bg, format='%Y-%m-%d %H:%M:%S')

# Print parsed configuration for verification
print(f"Simulation start date: {date_simul_bg}")
print(f"Station: {Station}")
print(f"Parameter structure: {param_struct}")
print(f"Sensors: {sensors}")
print(f"Sigma: {sigma}")
print(f"Zones to invert: {zones_to_invert}")
print(f"Parameters to invert: {parameters_to_invert}")
print(f"Simulations to do: {simul_todo}")
print(f"Simulation path: {main_dir}")
print(f"Ginette directory: {dir_ginette}")

# Load input data files
param_table = pd.read_csv("param_table.txt", sep=",", header=0)  # Parameter table
obs_temp = pd.read_csv("Obs_temp_PT100_t.dat", sep=",", header=0)  # Observed temperature data

# Load MSE results from sensitivity analysis
output_path = "SENSI_" + Station + "/S_mse_simul.dat"
mse_table = pd.read_csv(output_path, sep=",", header=0)

# Generate joint posterior plots for each zone
for zone in zones_to_invert:
    # Create parameter column names for this zone
    param_cols = [f"{param}{zone}" for param in parameters_to_invert]
    filename = f"{Station}_zone{zone}.png"
    # Plot joint posterior distribution
    plot_joint_posterior(df=mse_table, param_cols=param_cols, cmap="viridis", filename=filename)


for param in parameters_to_invert:
    # Create parameter column names for this zone
    param_cols = [f"{param}{zone}" for zone in zones_to_invert]
    filename = f"{Station}_param{param}.png"
    # Plot joint posterior distribution
    plot_joint_posterior(df=mse_table, param_cols=param_cols, cmap="viridis", filename=filename)

# Generate sensor-specific posterior plots for each zone
for zone in zones_to_invert:
    param_cols = [f"{param}{zone}" for param in parameters_to_invert]
    for sensor in sensors:
        post_col = f"{sensor}_mse_posterior"  # Posterior column name
        filename = f"{Station}_zone{zone}_{sensor}_mse_posterior.png"  # Output filename
        # Plot joint posterior by sensor
        plot_joint_posterior_by_sensor(
            df=mse_table, 
            param_cols=param_cols, 
            posterior_col=post_col, 
            cmap="viridis", 
            filename=filename, 
            plot_title=f"Posterior joint - Zone {zone} - {sensor}"
        )


# Generate sensor-specific posterior plots for each zone
for param in parameters_to_invert:
    param_cols = [f"{param}{zone}" for zone in zones_to_invert]
    for sensor in sensors:
        post_col = f"{sensor}_mse_posterior"  # Posterior column name
        filename = f"{Station}_param{param}_{sensor}_mse_posterior.png"  # Output filename
        # Plot joint posterior by sensor
        plot_joint_posterior_by_sensor(
            df=mse_table, 
            param_cols=param_cols, 
            posterior_col=post_col, 
            cmap="viridis", 
            filename=filename, 
            plot_title=f"Posterior joint - Zone {zone} - {sensor}"
        )