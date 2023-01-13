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

# name HZ temperature
temp_file=paste0(temp_sensor,"_",point_name,"_",dg_day,"_",dg_month,"_",dg_year,".csv")
temp_data=fread(temp_file,header = T)
temp_data=temp_data[, Filter(function(x) any(!is.na(x)), .SD)]

if (ncol(temp_data)==6) {
  colnames(temp_data)=c('n','dates','T1','T2','T3','T4')
}
if (ncol(temp_data)==5) {
  colnames(temp_data)=c('n','dates','T2','T3','T4')
}
if (ncol(temp_data)==4) {
  colnames(temp_data)=c('n','dates','T3','T4')
}

# number obs PT100 in the hyporheic zone, the last (deeper) one is used as boundary condition
ntemp= ncol(temp_data)-2
tempHobbo=data.frame(temp_data)
tempHobbo <- tempHobbo[, colSums(is.na(tempHobbo)) == 0]

tempHobbo$dates = as.POSIXct(temp_data$dates,'%d/%m/%Y %H:%M:%S', tz = 'GMT')

# read data pressure and river temperature
river_file=paste0(pres_sensor,"_",point_name,"_",dg_day,"_",dg_month,"_",dg_year,".csv")
riverHobbo=fread(river_file,header = T)
riverHobbo=riverHobbo[, Filter(function(x) any(!is.na(x)), .SD)]
riverHobbo <- riverHobbo[,1:4]
colnames(riverHobbo)=c('n','dates','pressure_differential_m','temperature_stream_C')
riverHobbo$dates = as.POSIXct(riverHobbo$dates,'%d/%m/%Y %H:%M', tz = 'GMT')
nriver= ncol(riverHobbo)-2
presDiff=data.frame(riverHobbo$dates,riverHobbo$pressure_differential_m)
stream_temp=data.frame(riverHobbo$dates,riverHobbo$temperature_stream_C)

#---- add consider timeseries from date of initial conditions ----
ini_pres= as.POSIXct(presDiff[1,1],'%d/%m/%Y %H:%M',tz='GMT')
end_pres=as.POSIXct(presDiff[dim(presDiff)[1],1],'%d/%m/%Y %H:%M',tz='GMT')

timeInitial = as.POSIXct(tempHobbo$dates[1],'%d/%m/%Y %H:%M',tz='GMT')
timeFinal = as.POSIXct(tempHobbo$dates[dim(tempHobbo)[1]],'%d/%m/%Y %H:%M',tz='GMT')

### DEPEND DE cal_time !!!!
end_date = min(end_pres,timeFinal,max(timeInitial,ini_pres,ini_date_cal) + cal_time)
begin_date=max(timeInitial,ini_pres,ini_date_cal)
tempHobbo<-subset(tempHobbo, dates > begin_date )
tempHobbo<-subset(tempHobbo, dates < end_date )
presDiff<-subset(presDiff,riverHobbo.dates<end_date)
stream_temp<-subset(stream_temp,riverHobbo.dates<end_date)
presDiff<-subset(presDiff,riverHobbo.dates>begin_date)
stream_temp<-subset(stream_temp,riverHobbo.dates>begin_date)
# reference times and dates starting from max(timeInitial,ini_obs)
t_time = seq(from = begin_date,
             to = min(end_pres,timeFinal, end_date),
             by = as.difftime(Deltat, units = "secs"))

t_dates = seq(from = strptime(paste0(as.character(t_time[1],format="%Y-%m-%d")," 00:00"),format = "%Y-%m-%d %H:%M"),
              to = strptime(paste0(as.character(t_time[length(t_time)],format="%Y-%m-%d")," 00:00"),format = "%Y-%m-%d %H:%M"),
              by = as.difftime(1, units = "days"))

# points where interpolation is to take place
xOut = as.numeric(difftime(t_time,t_time[1],units = "secs"))

xInit1 = as.numeric(difftime(presDiff[,1], t_time[1],units = "secs"))
presDiffInterp = spline(x=xInit1, y=presDiff[,2], xout = xOut)$y
#plot(xOut,presDiffInterp,type='l',xlim=c(50000,855000))

xInit2 = as.numeric(difftime(stream_temp[,1], t_time[1], units = "secs"))
tempStreamInterp = spline(x=xInit2,y=stream_temp[,2],xout = xOut)$y
#plot(xOut,tempStreamInterp,type='l',xlim=c(50000,55000))

xInit3 = as.numeric(difftime(tempHobbo$dates, t_time[1], units = "secs"))
tempHobboInterp = array(0, dim = c(length(xOut), ntemp))

for (i in 1:ntemp){
  tempHobboInterp[,i] = spline(x = xInit3, y = tempHobbo[, i+2], xout = xOut)$y
}

#plot(xOut,tempHobboInterp[,1],type='l',xlim=c(50000,70000))
# points(xInit,tempHobbo[,1])

#----define data for Ginette model----
sensorDepths = vector(length = ntemp)
for (i in 1:ntemp){
  if(i==1) {
    sensorDepths[i]=0-as.integer(depth_PT100[i,])*Deltaz
  } else {
    sensorDepths[i]=0+sensorDepths[i-1]-as.integer(depth_PT100[i,])*Deltaz
  }
}

z = seq(from = -Deltaz/2,
        to = sensorDepths[length(sensorDepths)] + Deltaz/2,
        by=-Deltaz)

### initial conditions

#IC pressure
# (Ginette needs hydraulic head  [m]) dp=h[HZ]-h[riv]
presIC = approx(x = c(0,sensorDepths[length(sensorDepths)]),
                y = c(-presDiffInterp[1]*9810,(0-sensorDepths[length(sensorDepths)])*9810),xout = z)$y

#IC temperature
tempIC = approx(x = c(0,sensorDepths),
                y = c(tempStreamInterp[1],tempHobboInterp[1,]),
                xout = z)$y

### boundary conditions
presBC = cbind(-presDiffInterp,0)  # Ginette needs hydraulic head [m]
tempBC = cbind(tempStreamInterp,tempHobboInterp[,ntemp])

### temperature timeseries for inversion
tsInv = tempHobboInterp[,1:nPT100] 

#depths of the temperature timeseries for inversion
depthsInv=-sensorDepths[1:(length(sensorDepths)-1)]/Deltaz #en cm
#Delete old files
file_path_obs<-"./GINETTE_SENSI/OBS/"
f_obs<-list.files(file_path_obs)
file.remove(f_obs)
#---- save Data ----
tOut = as.numeric(difftime(t_time,t_time[1],units = "secs"))
for (i in 1:nPT100) {
  data=data.frame(tOut,tempHobboInterp[,i])
  write.table(data,file = paste0("./GINETTE_SENSI/OBS/Obs_temperature_maille",i,"_t.dat"),col.names = FALSE,row.names = FALSE)
}
write.table(presBC,file = "./GINETTE_SENSI/E_charge_t.dat", col.names = FALSE,row.names = FALSE)
write.table(tempBC,file = "./GINETTE_SENSI/E_temp_t.dat",col.names = FALSE,row.names = FALSE)
write.table(presIC,file = "./GINETTE_SENSI/E_pression_initiale.dat", col.names = FALSE,row.names = FALSE)
write.table(tempIC,file = "./GINETTE_SENSI/E_temperature_initiale.dat",col.names = FALSE,row.names = FALSE)
write.table(tOut[length(tOut)]/86400,"./GINETTE_SENSI/nitt.dat",col.names = FALSE,row.names = FALSE)
model=data.frame(sensorDepths[length(sensorDepths)]/Deltaz*-1,sensorDepths[ntemp]*-1,t(depthsInv[seq(1,length(depthsInv))]))
write.table(model,"./GINETTE_SENSI/model.dat",col.names = FALSE,row.names = FALSE)


# #
# #
# ####### Plot all the observed data to check if everything's correct
# #
# #
# 
# # Temperature data
# Temp <- cbind(xOut, tempStreamInterp, tempHobboInterp)
# 
# # Nommer colonnes en fonction du nombre de PT100
# T_colomn_names <- c("time", "Stream")
# 
# for (i in seq_len(ntemp)) {
#   T_colomn_names[i+2] <- paste0(sensorDepths[i], " m")
# }
# colnames(Temp) <- T_colomn_names
# 
# Temp <- as.data.frame(Temp)
# Temp$time <- Temp$time + ini_date
# 
# # Mise en forme du data frame contenant les temperatures avec la fonction melt avant de tracer
# Melted_obs_t <- reshape2::melt(Temp, id.var = "time")
# colnames(Melted_obs_t) = c("time", "Depth", "Temperature")
# 
# # Recupération des éléments nécessaires pour nommer correctement les graphes
# T_titre <- paste0("T_Cal_Check_", as.character(point_name))
# 
# # Define colors
# colpal <-  brewer.pal(6, "Dark2")
# 
# # Plot
# g_temp_ts <-
#   ggplot() +
#   geom_line(data = Melted_obs_t,
#             mapping = aes(x = time, y = Temperature, color = Depth))  +
#   scale_color_manual(values = c(colpal)) +
#   labs(x = "", y = "T (C)", color = "Depth", title = T_titre) +
#   scale_x_datetime(date_labels = " %d %b %y") +
#   theme_bw()
# 
# #Save plot
# # png(paste0("PLOT/Data_check/", T_titre, ".png"))
# # g_temp_ts
# # dev.off()
# 
# 
# # IL FAUT UN TRAITEMENT AVANT
# 
# ## Head differential data
# Press <- as.data.frame(cbind(xOut, presDiffInterp))
# 
# # Nommer colonnes
# P_colomn_names <- c("p_time", "Head_differential")
# colnames(Press) <- P_colomn_names
# 
# #adaptation date
# Press$p_time <- Press$p_time + ini_dateme <- Press$p_time + ini_date
# 
# 
# # Recupération des éléments nécessaires pour nommer correctement les graphes
# P_titre <- paste0("P_Cal_Check_", as.character(point_name))
# 
# #plot
# g_meas_head_differential <-
#   ggplot() +
#   geom_line(data = Press,
#             mapping = aes(x = p_time,y = Head_differential)) +
#   expand_limits(y = 0) +
#   geom_hline(mapping = aes(yintercept = 0),linetype = "dashed") +
#   labs(x="",y = expression(Delta*'H = H'['HZ'] *'- H'['riv'] * ' (in m)'), title = P_titre) +
#   scale_x_datetime(date_labels ="%d %b %y") +
#   theme_bw()
# 
# # #Save plot
# # png(paste0("PLOT/Data_check/", P_titre, ".png"))
# # g_meas_head_differential
# # dev.off()