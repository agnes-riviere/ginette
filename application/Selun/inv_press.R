#LIBRARY
#---- Interpolate timeseries to match the time discretization ----
library(stats) #for spline interpolation
library(lubridate)
library(stringr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(data.table)
library(readr)
setwd(dir = ".")
getwd()
path_mini_lomos = '/home/ariviere/Programmes/ginette/application/Selun/'
setwd(path_mini_lomos)

# read comm
#File_com = list.files(pattern = "test_")
File_com = "inversion.COMM"
Inversion_PT100 = "inversion_PT100.COMM"
depth_PT100 = read.csv(Inversion_PT100, sep = " ", header = FALSE) #écart entre les PT100 en cm

com = read.csv(File_com, sep = " ", header = FALSE)
point_name=com$V1
temp_sensor=com$V5
pres_sensor=com$V6
#initial date
dg_year =com$V2
ini_year = as.numeric(paste0('20',dg_year))
ini_month = com$V3
dg_month = sapply(ini_month, function(x) paste(paste(rep(0, 2 - nchar(x)), collapse = ""), x, sep = ""))
ini_day = com$V4
dg_day=sapply(ini_day, function(x) paste(paste(rep(0, 2 - nchar(x)), collapse = ""), x, sep = ""))
ini_date = paste0(ini_day, '/', ini_month, '/', ini_year)
ini_date = as.POSIXct(ini_date, '%d/%m/%Y', tz = 'GMT')
cal_time=com$V11

dg_year_cal =com$V12
ini_year_cal = as.numeric(paste0('20',dg_year_cal))
ini_month_cal = com$V13
dg_month_cal = sapply(ini_month_cal, function(x) paste(paste(rep(0, 2 - nchar(x)), collapse = ""), x, sep = ""))
ini_day_cal = com$V14
dg_day_cal=sapply(ini_day_cal, function(x) paste(paste(rep(0, 2 - nchar(x)), collapse = ""), x, sep = ""))
ini_date_cal = paste0(ini_day_cal, '/', ini_month_cal, '/', ini_year_cal)
ini_date_cal = as.POSIXct(ini_date_cal, '%d/%m/%Y', tz = 'GMT')


# #---- discretisation parameters for Ginette ----
Deltaz = 0.01 # [m]
Deltat = com$V7 # [s]

# number obs PT100 in the hyporheic zone to inv
nPT100 = com$V8

# Liste tous les fichiers dans le répertoire
fichiers <- list.files('./')

# read data pressure and river temperature
# La chaîne de caractères par laquelle le nom du fichier doit commencer
chaine_river_file <- paste0(pres_sensor,"_",point_name)


# Filtre les fichiers dont le nom commence par la chaîne de caractères
river_file <- grep(paste0("^", chaine_river_file), fichiers, value = TRUE)
#river_file=paste0(pres_sensor,"_",point_name,"_",dg_day,"_",dg_month,"_",dg_year,".csv")
riverHobbo=fread(river_file,header = T)
riverHobbo=riverHobbo[, Filter(function(x) any(!is.na(x)), .SD)]
riverHobbo <- riverHobbo[,1:4]
colnames(riverHobbo)=c('n','dates','pressure_differential_m','temperature_stream_C')
#riverHobbo$dates = as.POSIXct(riverHobbo$dates,'%d/%m/%Y %H:%M', tz = 'GMT')
nriver= ncol(riverHobbo)-2
riverHobbo$pressure_differential_m=riverHobbo$pressure_differential_m*-1
# Formatez la date dans le format spécifié
#riverHobbo$dates <- format(riverHobbo$dates, format = "%d/%m/%Y %H:%M")


# Spécifiez le nom du fichier CSV dans lequel vous souhaitez sauvegarder la colonne
nom_fichier_csv <- paste0('inv',chaine_river_file,'.csv')

# Sauvegardez la colonne dans le fichier CSV
write.csv(riverHobbo, file = nom_fichier_csv, row.names = FALSE,quote = F)

