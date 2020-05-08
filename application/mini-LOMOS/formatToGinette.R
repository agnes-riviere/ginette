#LIBRARY
#---- Interpolate timeseries to match the time discretization ----
library(stats) #for spline interpolation
library(lubridate)
library(stringr)
library(ggplot2)
library(stringr)
library(RColorBrewer)

#wd=paste0('/home/ariviere/Documents/Bassin-Orgeval/Donnee_Orgeval_Mines/processed_data_KC/HZ/',
#          namePoint)
#setwd(wd)
# read comm
File_com = list.files(pattern = "test_")
Inversion_PT100 = "inversion_PT100.COMM"
depth_PT100 = read.csv(Inversion_PT100, sep = " ", header = FALSE) #écart entre les PT100 en cm

#ne marchera pas si length(File_com) != 1
com = read.csv(File_com, sep = " ", header = FALSE)

#initial date
ini_year = str_sub(com[1, 1], 7, 9)
ini_year = as.numeric(paste0('20', ini_year))
ini_month = str_sub(com[1, 1], 4, 5)
ini_day = str_sub(com[1, 1], 1, 2)
ini_date = paste0(ini_day, '/', ini_month, '/', ini_year)
ini_date = as.POSIXct(ini_date, '%d/%m/%Y', tz = 'GMT')

# #---- discretisation parameters for Ginette ----
Deltaz = 0.01 # [m]
Deltat = as.integer(as.numeric_version(com[3, 1])) # [s]

# number PT100 in the hyporheic zone
nPT100 = as.integer(as.numeric_version(com[4, 1]))

# name HZ temperature
for (i in 1:nPT100)  {
  a = read.csv(as.character(com[4 + i, 1]), sep = " ", header = FALSE)
  b = as.POSIXct(paste(a[, 1], a[, 2]), '%d/%m/%Y %H:%M:%S', tz = 'GMT')
  dat_tmp = data.frame(b, a[, 3])
  nam <- str_remove(as.character(com[4 + i, 1]), ".dat")
  if (i == 1) {temp1 = dat_tmp}
  assign(nam, dat_tmp)
}

tempHobbo = data.frame(Date = temp1[, 1])

for (i in 1:nPT100)  {
  a = read.csv(as.character(com[4 + i, 1]), sep = " ", header = FALSE)
  tempHobbo[, i + 1] = a[, 3]
}

# number river data
nriver= as.integer(as.numeric_version(com[nPT100+4+1,1]))

# read data pressure and river temperature
for (i in 1:nriver)  {
  a = read.csv(as.character(com[5 + nPT100 + i, 1]), sep = " ", header = FALSE)
  b = as.POSIXct(paste(a[, 1], a[, 2]), '%d/%m/%Y %H:%M:%S', tz = 'GMT')
  dat_tmp = data.frame(b, a[, 3])
  name <- str_remove(as.character(com[5 + nPT100 + i, 1]), ".dat")
  assign(nam, dat_tmp)
  if (i==1) {presDiff = dat_tmp}
  if (i==2) {stream_temp = dat_tmp}
}

#---- add consider timeseries from date of initial conditions ----
ini_pres= as.POSIXct(presDiff[1,1],'%d/%m/%Y %H:%M',tz='GMT')
end_pres=as.POSIXct(presDiff[dim(presDiff)[1],1],'%d/%m/%Y %H:%M',tz='GMT')

timeInitial = as.POSIXct(temp1[1,1],'%d/%m/%Y %H:%M',tz='GMT')
timeFinal = as.POSIXct(temp1[dim(temp1)[1],1],'%d/%m/%Y %H:%M',tz='GMT')


# reference times and dates starting from max(timeInitial,ini_obs)
t_time = seq(from = max(timeInitial,ini_pres),
             to = min(end_pres,timeFinal),
             by = as.difftime(Deltat, units = "secs"))

t_dates = seq(from = strptime(paste0(as.character(t_time[1],format="%Y-%m-%d")," 00:00"),format = "%Y-%m-%d %H:%M"),
              to = strptime(paste0(as.character(t_time[length(t_time)],format="%Y-%m-%d")," 00:00"),format = "%Y-%m-%d %H:%M"),
              by = as.difftime(1, units = "days"))

# points where interpolation is to take place
xOut = as.numeric(difftime(t_time,t_time[1],units = "secs"))

xInit = as.numeric(difftime(presDiff[,1],t_time[1],units = "secs"))
presDiffInterp = spline(x=xInit,y=presDiff[,2],xout = xOut)$y
#plot(xOut,presDiffInterp,type='l',xlim=c(50000,855000),ylim=c(-0.4,0.05))

xInit = as.numeric(difftime(stream_temp[,1],t_time[1],units = "secs"))
tempStreamInterp = spline(x=xInit,y=stream_temp[,2],xout = xOut)$y
#plot(xOut,tempStreamInterp,type='l',xlim=c(50000,55000))

xInit = as.numeric(difftime(tempHobbo$Date,t_time[1],units = "secs"))
tempHobboInterp = array(0,dim=c(length(xOut),nPT100))
for (i in 1:nPT100){
  tempHobboInterp[,i] = spline(x=xInit,y=tempHobbo[,i+1],xout = xOut)$y
}
#plot(xOut,tempHobboInterp[,1],type='l',xlim=c(50000,70000))
# points(xInit,tempHobbo[,1])

#----define data for Ginette model----
sensorDepths = vector(length = nPT100)
for (i in 1:nPT100){
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
tempBC = cbind(tempStreamInterp,tempHobboInterp[,nPT100])

### temperature timeseries for inversion
tsInv = tempHobboInterp[,1:3] 

#depths of the temperature timeseries for inversion
depthsInv=-sensorDepths[1:(length(sensorDepths)-1)]/Deltaz #en cm

#---- save Data ----
tOut = as.numeric(difftime(t_time,t_time[1],units = "secs"))
for (i in 1:nPT100) {
  data=data.frame(tOut,tempHobboInterp[,i])
  write.table(data,file = paste0("Obs_temperature_maille",i,"_t.dat"),col.names = FALSE,row.names = FALSE)
}
write.table(presBC,file = "E_charge_t.dat", col.names = FALSE,row.names = FALSE)
write.table(tempBC,file = "E_temp_t.dat",col.names = FALSE,row.names = FALSE)
write.table(presIC,file = "E_pression_initiale.dat", col.names = FALSE,row.names = FALSE)
write.table(tempIC,file = "E_temperature_initiale.dat",col.names = FALSE,row.names = FALSE)
write.table(tOut[length(tOut)],"nitt.dat",col.names = FALSE,row.names = FALSE)



# Plot all the observed data to check if everything's correct
# Temperature data
Temp <- cbind(xOut, tempStreamInterp, tempHobboInterp)

# Nommer colonnes en fonction du nombre de PT100
T_colomn_names <- c("time", "Stream")

if (nPT100 >= 2) {
  for (i in seq(2, nPT100)) {
    depth_PT100[i, 1] = depth_PT100[i - 1, 1] + depth_PT100[i, 1]
  }

}
depth_PT100 = -depth_PT100 / 100

for (i in seq_len(nPT100)) {
  T_colomn_names[i+2] <- paste0(depth_PT100[i, 1], " m")
}
colnames(Temp) <- T_colomn_names

# Mise en forme du data frame contenant les temperatures avec la fonction melt avant de tracer
Melted_obs_t <- reshape2::melt(as.data.frame(Temp), id.var = "time")
colnames(Melted_obs_t) = c("time", "Depth", "Temperature")
Melted_obs_t$time = Melted_obs_t$time + ini_date

# Recupération des éléments nécessaires pour nommer correctement les graphes
titre <- str_remove(File_com[1], pattern = ".COMM")
T_titre <- paste0("T_check_", str_remove(titre, pattern = "test_"))

# Define colors
colpal <-  brewer.pal(3, "Dark2")

# Plot
g_temp_ts <-
  ggplot() +
  geom_line(data = Melted_obs_t,
            mapping = aes(x = time, y = Temperature, color = Depth))  +
  scale_color_manual(values = c(colpal, colpal)) +
  labs(x = "", y = "T (C)", color = "Depth", title = "Temperature calibration Point 43") +
  scale_x_datetime(date_labels = " %d %b %y") +
  theme_bw()

#Save plot
png(paste0("PLOT/Data_check/", T_titre, ".png"))
g_temp_ts
dev.off()



## Head differential data
Press <- as.data.frame(cbind(xOut, presDiffInterp))

# Nommer colonnes
P_colomn_names <- c("p_time", "Head_differential")
colnames(Press) <- P_colomn_names

#adaptation date
Press$p_time <- Press$p_time + ini_date

# Recupération des éléments nécessaires pour nommer correctement les graphes
P_titre <- paste0("P_check_", str_remove(titre, pattern = "test_"))

#plot
g_meas_head_differential <-
  ggplot() +
  geom_line(data = Press,
            mapping = aes(x = p_time,y = Head_differential)) +
  expand_limits(y = 0) +
  geom_hline(mapping = aes(yintercept = 0),linetype = "dashed") +
  labs(x="",y = expression(Delta*'H = H'['HZ'] *'- H'['riv'] * ' (in m)'), title = "Head differential calibration Point43") +
  scale_x_datetime(date_labels ="%d %b %y") +
  theme_bw()

#Save plot
png(paste0("PLOT/Data_check/", P_titre, ".png"))
g_meas_head_differential
dev.off()