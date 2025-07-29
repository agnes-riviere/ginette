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


# Load MSE results from sensitivity analysis
output_path = "SENSI_" + Station + "/S_mse_simul.dat"

print("Looking for:", os.path.abspath(output_path))
print("Exists?", os.path.exists(output_path))

# print the path
print('path sensi',output_path)
mse_table = pd.read_csv(output_path, sep=",", header=0)
# Define parameter units (adjust according to your parameters)
param_units = {
    'k': 'm²',        # intrinsic permeability (square meters)
    'n': '-',         # porosity (dimensionless)
    'l': 'W/(m·K)',   # thermal conductivity (watts per meter-kelvin)
    # Add other parameters as needed
}

# Grid search analysis by zone
# This section creates parameter space plots for each geological zone being inverted
for zone in zones_to_invert:
    print(f"Grid Search for zone {zone}...")
    
    # Create parameter column names for this zone
    # Format: parameter_name + zone_number (e.g., "k1", "n1", "l1" for zone 1)
    param_cols = [f"{param}{zone}" for param in parameters_to_invert]
    filename = f"{Station}_zone{zone}_grid.png"
    
    # Plot grid search results
    # This creates a joint posterior plot showing parameter correlations and optimal values
    # - Shows 2D parameter space with color-coded MSE values
    # - Darker colors (viridis_r) indicate better model fits (lower MSE)
    # - Marginal distributions on axes show parameter sensitivity
    best_model = plot_joint_posterior(
        df=mse_table, 
        param_cols=param_cols, 
        filename=filename,
        param_units=param_units,
        cmap="viridis_r"  # Reversed viridis: dark = low MSE (good fit)
    )
    print(f"Grid search saved: {filename}")

print("Best model parameters:", best_model)

# Optional: Grid search by parameter
# This section creates plots showing how each parameter varies across all zones
for param in parameters_to_invert:
    print(f"Grid Search for parameter {param}...")
    
    # Create column names for this parameter across all zones
    # Shows spatial variability of each parameter (e.g., k1, k2, k3 for permeability)
    param_cols = [f"{param}{zone}" for zone in zones_to_invert]
    filename = f"{Station}_param{param}_grid.png"
    
    # Plot parameter correlations across zones
    # Helps identify if parameters are spatially correlated or independent
    # - Strong correlations suggest similar hydrogeological properties
    # - Weak correlations indicate heterogeneous subsurface conditions
    best_model = plot_joint_posterior(
        df=mse_table, 
        param_cols=param_cols, 
        filename=filename,
        param_units=param_units,
        cmap="viridis_r"
    )
    print(f"Grid search saved: {filename}")

print("Best model parameters:", best_model)



# Generate sensor-specific posterior plots for each zone
# This analysis shows how well each sensor is fitted by different parameter combinations
for zone in zones_to_invert:
    param_cols = [f"{param}{zone}" for param in parameters_to_invert]
    
    for sensor in sensors:
        post_col = f"{sensor}"  # MSE column for this specific sensor
        filename = f"{Station}_zone{zone}_{sensor}_mse.png"
        
        # Plot parameter space colored by individual sensor MSE
        # Interpretation guide:
        # - Each sensor responds differently to parameter changes
        # - Temperature sensors: sensitive to thermal conductivity and porosity
        # - Pressure sensors: sensitive to permeability and storage coefficients
        # - Multi-objective optimization balances all sensor responses
        plot_joint_posterior_by_sensor(
            df=mse_table, 
            param_cols=param_cols, 
            col=post_col, 
            cmap="viridis_r", 
            filename=filename
        )


# Generate sensor-specific posterior plots for each zone
for param in parameters_to_invert:
    param_cols = [f"{param}{zone}" for zone in zones_to_invert]
    for sensor in sensors:
        post_col = f"{sensor}"  # Posterior column name
        filename = f"{Station}_param{param}_{sensor}_mse.png"  # Output filename
        # Plot joint posterior by sensor
        plot_joint_posterior_by_sensor(
            df=mse_table, 
            param_cols=param_cols, 
            col=post_col, 
            cmap="viridis_r", 
            filename=filename
        )
