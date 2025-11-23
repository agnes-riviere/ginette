#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 14:00:41 2025

@author: Maxime Gautier
"""


# %% IMPORTS:
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt


# %% FUNCTIONS:
def misfit(obs, sim, err):
    return np.sum(((sim - obs) / err)**2)


def likelyhood(m):
    return np.exp(-0.5*m)


def log_likelyhood(m):
    return np.log(likelyhood(m))

# %% MISFIT:
results = pd.read_csv(os.path.join("grid_search.csv"), delimiter=";")
obs_data = pd.read_csv(os.path.join("observed_data.txt"), delimiter=" ",
                       index_col=[0])
results["misfit_1"] = np.nan
results["misfit_2"] = np.nan
results["misfit_3"] = np.nan

done = sorted([int(f.split("_")[-1].split(".")[0])
               for f in os.listdir(os.path.join("results"))])[1:]

# Misfit by simulations:
err = 0.2
results["err_misfit"] = err
for i in tqdm(done, desc="Compute misfit"):
    sim = pd.read_csv(os.path.join("results", "", f"sim_temp_{i}.txt"),
                      delimiter=" ", index_col=[0])
    m1 = misfit(obs_data.Temp1, sim.Temp1, err=err)
    m2 = misfit(obs_data.Temp2, sim.Temp2, err=err)
    m3 = misfit(obs_data.Temp3, sim.Temp3, err=err)

    results.loc[results["ID"] == i, "misfit_1"] = m1
    results.loc[results["ID"] == i, "misfit_2"] = m2
    results.loc[results["ID"] == i, "misfit_3"] = m3

# Total misfit and likelihood
results["misfit_tot"] = (results["misfit_1"] + results["misfit_2"] +
                         results["misfit_3"])

results["L1"] = likelyhood(results["misfit_1"])
results["L2"] = likelyhood(results["misfit_2"])
results["L3"] = likelyhood(results["misfit_3"])

results.to_csv(os.path.join("results.txt"), sep=" ")