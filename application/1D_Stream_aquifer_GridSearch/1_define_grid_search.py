#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 20:19:17 2024

@author: Maxime GAUTIER
"""


# %% IMPORTS:
import os
import itertools
import numpy as np
import pandas as pd


# %% PARAMETERS:
# Define ranges for test parameters in the grid search:
N = 100  # number of test by parameter (nb)

"""
# True values:
log_k = -12  # log of intrinsec permeability (k in m^2)
lam = 2.5  # thermal conductivity (W/m.°C)
"""

# Parameters ranges [min, max]:
log_k = [-15, -11]  # log of intrinsec permeability (k in m^2)
lam = [1, 8]  # thermal conductivity (W/m.°C)


# %% CREATE TABLE FOR THE GRID SEARCH:
log_k_grid = np.logspace(log_k[0], log_k[1], N)
lam_grid = np.linspace(lam[0], lam[1], N)

grid_search = [t for t in itertools.product(log_k_grid, lam_grid)]

grid_search = pd.DataFrame(grid_search, columns=["log_k", "lam"])
grid_search.index.name = "ID"


# %% SAVE:
grid_search.to_csv(os.path.join("grid_search.csv"), sep=";")
