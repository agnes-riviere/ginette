#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08-12-2025


@author: Agnès Rivière, Samuel Larance
"""


# IMPORT:
import importlib
import os
import sys
import numpy as np
import pandas as pd
from time import time
import shutil
import multiprocessing as mp
from pathlib import Path
import glob
from src.src_python import Read_obs



# Set up directories and data paths

# Get current script directory
SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = str(SCRIPT_DIR)

# Navigate to sibling directories
GINETTE_SENSI = os.path.join(BASE_DIR, "GINETTE_SENSI")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
# This path contains the observational data (temperature sensors, pressure measurements)
Obs_data = os.path.join(BASE_DIR,'OBS_point/Point3_540_SOULTZ/')
# =============================================================================
# SIMULATION SETUP: Define temporal parameters and read observational data
# =============================================================================
importlib.reload(Read_obs)
# Start date and time for the simulation
# This corresponds to when field measurements began
date_simul_bg = pd.to_datetime("2025/09/26 12:00:00")

# Simulation state configuration:
# 0 = steady state (time-independent, equilibrium conditions)
# 1 = transient state (time-dependent, dynamic evolution)
# We use transient state to capture temporal variations in temperature and flow
state = 1

# Data processing coefficients
coef = 1    # Scaling coefficient for pressure measurements
offset = 0  # Offset for pressure measurements (baseline correction)
# Simulation duration in days
# 30 days provides sufficient time to observe seasonal patterns and model spin-up
nb_day = 30
# =============================================================================
# MODEL GEOMETRY AND DISCRETIZATION SETUP
# =============================================================================

# TEMPORAL DISCRETIZATION
# Time step in seconds (900s = 15 minutes)
# This matches the measurement frequency and ensures numerical stability
dt = 900

# SPATIAL DOMAIN DEFINITION (1D vertical column)
z_top = 0.0      # Surface elevation (stream bed) [m]
z_bottom = -0.4  # Bottom of model domain [m] (40 cm below streambed)
az = abs(z_top - z_bottom)  # Total column height [m]

# GRID DISCRETIZATION
dz = 0.01        # Vertical cell size [m] (1 cm resolution)
                 # Fine discretization needed to capture thermal gradients

# OBSERVATION DEPTHS
dz_obs = 0.1     # Spacing between temperature sensors [m] (10 cm)
                 # Sensors at -10, -20, -30, -40 cm depths

# =============================================================================
# HYDROGEOLOGICAL PARAMETER DEFINITION
# =============================================================================

# GEOLOGICAL HETEROGENEITY
# Number of geological facies (material zones) in the column
# nb_zone = 1: Homogeneous porous medium
# nb_zone = 2: Two-layer system (typical for streambed environments)
nb_zone = 1

# Boundary between geological zones [m]
# Negative value indicates depth below streambed surface
alt_thk = -0.32  # Interface at 32 cm depth

# =============================================================================
# ZONE 1 PARAMETERS (Upper layer: 0 to -32 cm)
# Typically represents looser, more permeable streambed sediments
# =============================================================================

# POROSITY (φ) - Fraction of void space available for fluid flow
REF_n = 0.8  # High porosity typical of loose gravel/sand

# INTRINSIC PERMEABILITY (k) - Measure of medium's ability to transmit fluid
# Relationship: k = K·μ/(ρ·g) where:
# - K: hydraulic conductivity [m/s]
# - μ: dynamic viscosity [Pa·s]
# - ρ: fluid density [kg/m³]
# - g: gravitational acceleration [m/s²]
REF_k = -15  # Log₁₀(permeability in m²)
               # k = 10^(-13.5) ≈ 3.16×10^(-14) m²
               # Corresponds to coarse sand/fine gravel

# THERMAL CONDUCTIVITY (λ) - Heat conduction efficiency [W/m·K]
REF_l = 2.0    # Typical for water-saturated sediments

# SOLID GRAIN DENSITY (ρₛ) - Density of mineral grains [kg/m³]
REF_r = 3500   # Typical for quartz-rich sediments

# HEAT CAPACITY CALCULATION (handled internally by Ginette):
# c_pm = c_w·ρ_w·n·S + c_s·ρ_s·(1-n) + c_a·ρ_a·n·(1-S)
# where:
# - c_w = 4185 J/kg·°C (specific heat of water)
# - c_s = variable J/kg·°C (specific heat of solid)
# - S = saturation ratio
# - Fixed ρ_m = 1000 kg/m³ (imposed by Ginette)

# =============================================================================
# ZONE 2 PARAMETERS (Lower layer: -32 to -40 cm)
# Typically represents more compact, less permeable substrate
# =============================================================================

if nb_zone == 2:
    REF_n2 = 0.3    # Lower porosity (more compact sediment)
    REF_k2 = -13     # Higher permeability: k = 10^(-13) m²
    REF_l2 = 2.65    # Higher thermal conductivity (more mineral content)
    REF_r2 = 3500    # Same grain density (similar mineral composition)

print(f"Geological configuration:")
print(f"- Number of zones: {nb_zone}")
if nb_zone == 2:
    print(f"- Zone interface at: {alt_thk} m")
    print(f"- Zone 1 (upper): n={REF_n}, k=10^{REF_k} m², λ={REF_l} W/m·K")
    print(f"- Zone 2 (lower): n={REF_n2}, k=10^{REF_k2} m², λ={REF_l2} W/m·K")
else:
    print(f"- Homogeneous medium: n={REF_n}, k=10^{REF_k} m², λ={REF_l} W/m·K")




# Find project-relative src and application directories (no absolute paths)
def find_project_paths(start_file=__file__):
    p = Path(start_file).resolve().parent
    repo_root = None
    src_py = None
    app_dir = None
    for _ in range(8):
        # candidate locations relative to current parent
        cand_src = p / "src" / "src_python"
        cand_src2 = p / "src_python"
        cand_app = p / "application" / "1D_Stream_aquifer_GridSearch"
        if cand_src.exists():
            src_py = str(cand_src)
        if cand_src2.exists() and src_py is None:
            src_py = str(cand_src2)
        if cand_app.exists():
            app_dir = str(cand_app)
        if src_py or app_dir:
            repo_root = str(p)
            break
        p = p.parent
    return repo_root, src_py, app_dir

REPO_ROOT, SRC_PY, BASE_APP_DIR = find_project_paths()

# Fallback to reasonable relative defaults if not found
if REPO_ROOT is None:
    REPO_ROOT = Path(__file__).resolve().parents[2].as_posix()
if SRC_PY is None:
    SRC_PY = os.path.join(REPO_ROOT, "src", "src_python")
if BASE_APP_DIR is None:
    BASE_APP_DIR = os.path.join(REPO_ROOT, "application", "1D_Stream_aquifer_GridSearch")

# add src python to sys.path if available
if os.path.isdir(SRC_PY) and SRC_PY not in sys.path:
    sys.path.insert(0, SRC_PY)
if BASE_APP_DIR is None:
    BASE_APP_DIR = os.path.join(REPO_ROOT, "application", "1D_Stream_aquifer_GridSearch")

# --- Ajout: dossier de sortie local à l'application (modifiable) ---
RESULTS_DIR = os.path.join(BASE_APP_DIR, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# debug listing to verify (non-blocking)
if os.path.isdir(SRC_PY):
    try:
        print(f"Using src_python: {SRC_PY}")
        print("Files:", sorted(os.listdir(SRC_PY)))
    except Exception:
        pass

# Import project modules from src_python (robust to module name/case)
try:
    # preferred: modules as they appear in src/src_python
    from Direct_model import (setup_ginette2,setup_ginette,
                               initial_conditions,
                               boundary_conditions,
                               run_direct_model,
                               smooth_square_wave)
    from Read_obs import (process_obs_data, convert_dates,read_csv_with_multiple_separators)
    from Plot import plot_initial_conditions
except Exception:
     # fallback to legacy module name if present
    from Direct_model import (setup_ginette2,
                                       initial_conditions,
                                       boundary_conditions,
                                       run_direct_model,
                                       smooth_square_wave)


try:
    from Init_folders import prepare_ginette_directories
except Exception:
    # fallback if an alternative utils module exists
    try:
        from utils_ginette import prepare_ginette_directories
    except Exception:
        # minimal fallback: create directory helper
        def prepare_ginette_directories(path):
            os.makedirs(path, exist_ok=True)

# provide a portable copy_file helper if the project does not expose one
try:
    from Init_folders import copy_file  # some versions may provide it
except Exception:
    def copy_file(src, dst_dir):
        """Copy src into dst_dir (create dst_dir if needed)."""
        os.makedirs(dst_dir, exist_ok=True)
        if not os.path.exists(src):
            raise FileNotFoundError(f"Source not found: {src}")
        shutil.copy(src, os.path.join(dst_dir, os.path.basename(src)))







# List available CSV files in the observational data directory
print("CSV files in observational data directory:")
for file in glob.glob(os.path.join(Obs_data, '*.csv')):
    print(f"- {file}")

# Process observational data using the corrected function:
# - Reads temperature and pressure time series from CSV files with proper semicolon parsing
# - Interpolates to 15-minute intervals (matching model time step)
# - Applies quality control and filtering
# - Returns DataFrame with synchronized measurements
obs_temp = process_obs_data(Obs_data, date_simul_bg, coef, offset, nb_day)
# modifier pour le pas de simulation soit variable et pas 900 secondes
print(f"\nObservational data loaded successfully:")
print(f"- Time period: {obs_temp.index.min()} to {obs_temp.index.max()}")
print(f"- Number of time steps: {len(obs_temp)}")
print(f"- Available measurements: {list(obs_temp.columns)}")
print(f"- Data shape: {obs_temp.shape}")


# GO in  BASE_DIR
os.chdir(GINETTE_SENSI)
# Initialize Ginette model files and return observation depths
# This function creates all necessary input files for the Ginette model:
# - Parameter files (E_parametre.dat)
# - Thermal parameter files (E_p_therm.dat)
# - Grid coordinates and observation points
z_obs = setup_ginette(dt, state, nb_day, z_top, z_bottom, az, dz, date_simul_bg, dz_obs)

print(f"Model domain setup:")
print(f"- Vertical extent: {z_top} to {z_bottom} m")
print(f"- Cell size: {dz} m")
print(f"- Number of cells: {int(az/dz)}")
print(f"- Time step: {dt} s ({dt/60} minutes)")
print(f"- Observation depths: {z_obs} m")



# =============================================================================
# INITIAL AND BOUNDARY CONDITIONS SETUP
# =============================================================================

# INITIAL CONDITIONS: Set starting temperature and pressure fields
# These represent the system state at t=0 (simulation start)
initial_conditions(obs_temp, z_top, z_bottom, dz, z_obs)
# This function:
# 1. Interpolates measured temperatures to all grid cells
# 2. Creates initial pressure profile from surface pressure measurements
# 3. Writes E_temperature_initiale.dat and E_charge_initiale.dat

# BOUNDARY CONDITIONS: Define time-varying surface and bottom conditions
# Critical for realistic simulation of stream-aquifer interactions
boundary_conditions(obs_temp, dt)
# This function:
# 1. Sets surface temperature from stream measurements (time-varying)
# 2. Sets bottom temperature from deep sensor (less variable)
# 3. Sets surface pressure from differential pressure measurements
# 4. Sets bottom pressure boundary (typically fixed or gradient)
# 5. Writes E_cdt_aux_limites.dat, E_charge_t.dat, E_temp_t.dat

print("Initial and boundary conditions applied:")
print("- Initial temperature profile interpolated from observations")
print("- Initial pressure profile derived from surface measurements")
print("- Time-varying surface conditions from observational data")


# =============================================================================
# PLOT INITIAL CONDITIONS
# =============================================================================



# Plot the initial temperature and pressure profiles
# This visualization helps verify that the initial conditions are realistic
# and properly interpolated from the observational data
plot_initial_conditions()

# The plots show:
# 1. Initial temperature profile: Interpolated from sensor measurements
# 2. Initial pressure profile: Derived from surface pressure measurements
# 
# These initial conditions serve as the starting point (t=0) for the simulation
# and should represent a reasonable approximation of the actual field conditions
