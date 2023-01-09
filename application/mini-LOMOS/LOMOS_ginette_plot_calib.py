#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import subprocess
import rpy2.robjects as robjects

path_mini_lomos='/home/ariviere/Programmes/ginette/application/mini-LOMOS/'
os.chdir(path_mini_lomos)

path_plot='PLOT'
os.chdir(os.path.join(path_mini_lomos,path_plot))
import subprocess
namePointT='T3_Point3Nelly_14_04_22.csv'
namePointP='P2_Point3Nelly_14_04_22.csv'
date_begin='14/04/2022 17:45:00'
sim_name=1

#functions R# - 
robjects.r.source("Plot_function.R")

#Plot simulated and measured Temperature time series
robjects.r['temperature_ts'](sim_name,date_begin)

#Plot Temperature profile
#robjects.r['temperature_profile'](sim_name)



#Plot conductive_fluxes interpolated with time and depth
#robjects.r['conductive_fluxes'](sim_name)

#Plot advective_fluxes
#robjects.r['advective_fluxes'](sim_name)


#Plot fluxes interface stream-aquifer
#robjects.r['flux_ts'](namePointT,namePointP,sim_name,date_begin)


#All profile with date
#robjects.r['tot_fig'](namePointT,namePointP,sim_name,date_begin)

