# set working directory (auto-detected: works with Rscript, and with RStudio
# when the rstudioapi package is installed; falls back to the current
# working directory otherwise)
find_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd_args[grepl("^--file=", cmd_args)])
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(file_arg)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }
  getwd()
}
path_mini_lomos = paste0(find_script_dir(), "/")
setwd(path_mini_lomos)

# set plot directory
path_plot = 'PLOT'
setwd(paste(path_mini_lomos, path_plot, sep = ""))

# set input file names
namePointT = 't5_p5amont1_12_05_21.csv'
namePointP = 'p5_p5amont1_12_05_21.csv'

# set date
date_begin = '23/05/2021 17:00:00'

# set simulation name
sim_name = 1

# source R script with functions
source("Plot_function.R")

# plot temperature time series
temperature_ts(sim_name, date_begin)

# plot temperature profile
#temperature_profile(sim_name)

# plot conductive fluxes
#conductive_fluxes(sim_name)

# plot advective fluxes
#advective_fluxes(sim_name)

# plot fluxes at interface between stream and aquifer
flux_ts(namePointT, namePointP, sim_name, date_begin)

# plot all profiles with date
tot_fig(namePointT, namePointP, sim_name, date_begin)

