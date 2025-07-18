# -*- coding: utf-8 -*-
import os
import sys
# RIV2D + Station
dir_ginette = "/home/ariviere/Programmes/ginette"
Station = "AmB"  # Example station, replace with actual station name
sys.path.append(dir_ginette)  # Ajouter le dossier parent de src

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

# Ajouter le chemin du dossier src au PYTHONPATH
import importlib

from src.src_gmsh import mesh_generator
from src.src_python import Init_folders
from src.src_python import Direct_model
from src.src_python import Read_obs
from src.src_python import Plot
from src.src_python import stat_critere
from src.src_python import Grid_search

# Import all functions/classes from the relevant modules
from src.src_gmsh.mesh_generator import *
from src.src_python.Init_folders import *
from src.src_python.Direct_model import *
from src.src_python.Read_obs import *
from src.src_python.Plot import *
from src.src_python.stat_critere import *
from src.src_python.Grid_search import *



importlib.reload(Init_folders)
importlib.reload(Direct_model)
importlib.reload(Grid_search)
from src.src_python.Grid_search import *
from src.src_python import Grid_search
importlib.reload(Grid_search)


importlib.reload(Init_folders)
importlib.reload(Direct_model)
importlib.reload(Grid_search)

print('repertoire:',dir_ginette)

main_dir = os.path.join(dir_ginette, "application/RIV2D", Station)
os.chdir(main_dir)
print(main_dir)
with open("input_inversion.txt", "r") as f:
    lines = f.readlines()

# Helper to parse a line of the form "key: value"
def parse_line(line):
    key, value = line.strip().split(":", 1)
    return key.strip(), value.strip()

parsed = dict(parse_line(line) for line in lines)

date_simul_bg = parsed["date_simul_bg"]
Station = parsed["Station"]
# Read param_struct from the parsed input and convert to list
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




# Convert the date_simul_bg to datetime
date_simul_bg = pd.to_datetime(date_simul_bg, format='%Y-%m-%d %H:%M:%S')
print(f"date_simul_bg: {date_simul_bg}")
print(f"Station: {Station}")
print(f"param_struct: {param_struct}")
print(f"sensors: {sensors}")
print(f"sigma: {sigma}")
print(f"zones_to_invert: {zones_to_invert}")
print(f"parameters_to_invert: {parameters_to_invert}")
print(f"simul_todo: {simul_todo}")
print(f"path_simul: {main_dir}")
print(f"dir_ginette: {dir_ginette}")

# read the param_table.txt file
param_table = pd.read_csv("param_table.txt", sep=",", header=0)

# Load observed temperature data from CSV file
obs_temp = pd.read_csv("Obs_temp_PT100_t.dat", sep=",", header=0)

output_path = "SENSI_" + Station + "/S_mse_simul.dat"
mse_table = pd.read_csv(output_path, sep=",", header=0)


#print(f"DataFrame enregistré dans : {output_path}")

#Sortir meilleur model
best_row = mse_table.loc[mse_table['Total_mse'].idxmin()]
# build params_of_interest with parameters_to_invert  zones_to_invert 

params_of_interest = best_row[["k4", "n4", "l4", "k5", "n5", "l5", 'Total_mse']]
print("Valeurs des paramètres pour le plus petit Total_mse :")
print(params_of_interest)

#Sortir meilleur model
best_row = mse_table.loc[mse_table['Total_mse'].idxmin()]
params_of_interest = best_row[["k4", "n4", "l4", "k5", "n5", "l5", 'Total_mse']]
print("Valeurs des paramètres pour le plus petit Total_mse :")
print(params_of_interest)

for zone in zones_to_invert:
    param_cols = [f"{param}{zone}" for param in parameters_to_invert]
    filename = f"{Station}_posterior_zone{zone}.png"
    plot_joint_posterior(df=mse_table,param_cols=param_cols,cmap="viridis",filename=filename)


for zone in zones_to_invert:
    param_cols = [f"{param}{zone}" for param in parameters_to_invert]
    for sensor in sensors:
        post_col = f"{sensor}_mse_posterior" # <- colonne des postérieurs
        filename = f"{Station}_zone{zone}_{sensor}_mse_posterior.png" # <- nom du fichier
        plot_joint_posterior_by_sensor(df=mse_table,param_cols=param_cols,posterior_col=post_col,cmap="viridis",filename=filename,plot_title=f"Posterior joint - Zone {zone} - {sensor}")

