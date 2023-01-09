# set working directory
path_mini_lomos = '/home/ariviere/Programmes/ginette/application/mini-LOMOS/'
setwd(path_mini_lomos)

# set plot directory
path_plot = 'PLOT'
setwd(paste(path_mini_lomos, path_plot, sep = ""))

# set input file names
namePointT = 'T3_Point3Nelly_14_04_22.csv'
namePointP = 'P2_Point3Nelly_14_04_22.csv'

# set date
date_begin = '24/05/2022 17:45:00'

# set simulation name
sim_name = 2

# source R script with functions
source("Plot_function.R")

# plot temperature time series
temperature_ts(sim_name, date_begin)

# plot temperature profile
temperature_profile(sim_name)

# plot conductive fluxes
conductive_fluxes(sim_name)

# plot advective fluxes
advective_fluxes(sim_name)

# plot fluxes at interface between stream and aquifer
flux_ts(namePointT, namePointP, sim_name, date_begin)

# plot all profiles with date
tot_fig(namePointT, namePointP, sim_name, date_begin)

