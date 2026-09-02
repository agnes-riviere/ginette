library(data.table)
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


File_com = "inversion.COMM"
# set input file names
com = read.csv(File_com, sep = " ", header = FALSE)
# Obtenir le répertoire de travail
repertoire_de_travail <- getwd()

# Imprimer le répertoire de travail
print(repertoire_de_travail)


# set plot directory
path_plot = 'PLOT'
setwd(paste(path_mini_lomos, path_plot, sep = ""))


point_name=com$V1
temp_sensor=com$V5
pres_sensor=com$V6
# name HZ temperature

# La chaîne de caractères par laquelle le nom du fichier doit commencer
chaine_temp_file <- paste0(temp_sensor,"_",point_name)

# Liste tous les fichiers dans le répertoire
fichiers <- list.files('../')

# Filtre les fichiers dont le nom commence par la chaîne de caractères
namePointT <- grep(paste0("^", chaine_temp_file), fichiers, value = TRUE)


# read data pressure and river temperature
# La chaîne de caractères par laquelle le nom du fichier doit commencer
chaine_river_file <- paste0(pres_sensor,"_",point_name)

# Liste tous les fichiers dans le répertoire

# Filtre les fichiers dont le nom commence par la chaîne de caractères
namePointP <- grep(paste0("^", chaine_river_file), fichiers, value = TRUE)




print(namePointT)
print(namePointP)
# set date
#initial date
dg_year =com$V12
ini_year = as.numeric(paste0('20',dg_year))
ini_month = com$V13
dg_month = sapply(ini_month, function(x) paste(paste(rep(0, 2 - nchar(x)), collapse = ""), x, sep = ""))
ini_day = com$V14
dg_day=sapply(ini_day, function(x) paste(paste(rep(0, 2 - nchar(x)), collapse = ""), x, sep = ""))
ini_date = paste0(ini_day, '/', ini_month, '/', ini_year)
date_begin = as.POSIXct(ini_date, '%d/%m/%Y', tz = 'GMT')


# set simulation name
sim_name = 4

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

