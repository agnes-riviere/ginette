## plot velocity time series in cm/day
library(data.table)
namePointT <- 'T3_Point3Nelly_14_04_22.csv'
namePointP <- 'P2_Point3Nelly_14_04_22.csv'
iFile=58
sim_name=1

date_begin=as.POSIXct('14/04/2022 17:45:00',tz = 'GMT',format='%d/%m/%Y %H:%M:%S')


#paths
path_output <- "../GINETTE_SENSI/OUTPUT/"
path_plot <- "../PLOT/"
path_obs <- "../GINETTE_SENSI/OBS/"
path_BC <- "../GINETTE_SENSI/"


# -- read files with velocity measurements --

#pathFile = paste0(namePoint,'/path_sim/iFile',
#                  'OUTPUT/S_vitesse_nmaille2_hb',isim,'.dat')
pathFile_velocity =paste0(path_output,'S_vitesse_nmaille2_hb_',sim_name,'.dat')


data_velocity_raw <- fread(file = pathFile_velocity,sep=' ',header = F)
data_velocity <- data.frame(idx=data_velocity_raw[,1],
                            v1=data_velocity_raw[,2] * (100*3600*24), # from m/s to cm/d
                            v2=data_velocity_raw[,3]* (100*3600*24)) # from m/s to cm/d
rm(data_velocity_raw)

# -- read files with flux --

pathFile_heat_flx = paste0(path_output,'S_flux_therm_velocity_1_t_',sim_name,'.dat')

data_flux_raw_temp <- fread(file = pathFile_heat_flx,sep=' ',header = F)

# colonne 1 : date
# colonne 2 : conduction
# colonne 3 : advection
# colonne 4 : somme
# colonne 5 : darcy
data_flux_raw=data.frame()
data_flux_raw <- data.frame(data_flux_raw_temp[,1:5])
colnames(data_flux_raw) <- c('time','cond','adv','sum','darcy')

# changer la vitesse de Darcy de m/s a cm/day 
data_flux_raw$darcy <- data_flux_raw$darcy * (100*3600*24) # from m/s to cm/d
data_flux_raw$dates<-data_flux_raw$time+date_begin


# ---- plot all together ----

time1 = date_begin
time2 = max(data_flux_raw$dates,na.rm=T)
# -- read files to get dates --

pathMeas = paste0(path_BC,'../',namePointT)
data_meas <- fread(file = pathMeas,header = T,sep = ',')
colnames(data_meas)=c('id','dates','T_depth_10cm_C','T_depth_20cm_C','T_depth_30cm_C','T_depth_40cm_C')
data_meas$dates <- as.POSIXct(x = data_meas$dates,tz = 'GMT',format='%d/%m/%Y %H:%M:%S')
data_meas=data_meas[data_meas$dates >= time1&data_meas$dates <= time2, ]
pathMeasP = paste0(path_BC,'../',namePointP)
data_measP <- fread(file = pathMeasP,header = T,sep = ',')
colnames(data_measP)=c('id','dates','dh','T_stream_C')
data_measP$dates <- as.POSIXct(x = data_measP$dates,tz = 'GMT',format='%d/%m/%Y %H:%M:%S')

date_end=min(max(data_meas$dates),max(data_measP$dates))
dates_perDay <- seq(from=date_begin,
                    to=date_end,
                    by="day")
dates_perWeek <- dates_perDay[seq(from=1,to=length(dates_perDay),by=7)]
data_measP=data_measP[data_measP$dates >= time1 &data_measP$dates <= time2 , ]

## change later - 
# j'avais oublie de changer le nombre de pas de temps dans ginette
# a faire automatiquement
data_velocity <- data_velocity[1:nrow(data_meas),]
data_flux_raw <- data_flux_raw[1:nrow(data_meas),]





# in the plot we want to add
# whether the temperature gradient in the river
# is positive or negative

# to do so, check the temperature timeseries
# in the river and at the bottom of the HZ
data_measP$T_stream_C <- as.numeric(data_measP$T_stream_C)
data_meas$T_depth_40cm_C <- as.numeric(data_meas$T_depth_40cm_C)

plot(data_measP$T_stream_C,type='l')
lines(data_meas$T_depth_40cm_C)

# isPosTherm contains true if the thermal gradient is positive
# (ie warmer in the stream than in the hyporheic zone)
data_meas$isPosTherm <- ((data_measP$T_stream_C - data_meas$T_depth_40cm_C) > 0)
idx_pos2neg <- 1 + which(diff(data_meas$isPosTherm) == -1)
idx_neg2pos <- 1 + which(diff(data_meas$isPosTherm) == 1)
df_isPos <- data.frame(begin=c(idx_neg2pos), end=idx_pos2neg)

for(i in 1:nrow(df_isPos)){
  polygon(x=c(df_isPos[i,],rev(df_isPos[i,])),
          y=c(rep(0,2),rep(20,2)),
          col='peachpuff',border = NA
            )
}
lines(data_meas$T_stream_C)
lines(data_meas$T_depth_40cm_C)
dev.off()
# how to plot that now....
# 

png(file = paste0('../PLOT/flux_timeseries.png'),width = 1000,height = 500)
op<-par(oma=c(1,1,1,1),mar = c(1, 3, 0, 3)+0.1,cex=1,cex.main=1.5, cex.lab=1.5, cex.axis=1.2)

# plot Darcy velocity on first y-axis
plot(x = data_flux_raw$dates,y = data_flux_raw$darcy,type='l',
     xlab='',ylab='',xaxt='n',xlim=c(time1,time2),ylim=c(min(data_flux_raw$darcy,na.rm = T),max(data_flux_raw$darcy,na.rm = T)))
mtext(text = 'water velocity [cm/day]',side = 2,line = 2.5,cex=1.5)
axis.POSIXct(side = 1,x = data_flux_raw$dates,at=dates_perWeek,
             format = '%m/%d/%y')

# add background for positive thermal gradient
for(i in 1:nrow(df_isPos)){
  datesLim_i <- data_meas$dates[as.numeric(df_isPos[i,])]
  polygon(x=c(datesLim_i,rev(datesLim_i)),
          y=c(rep(-10,2),rep(10,2)),
          col='peachpuff',border = NA
  )
}

# add horizontal grid
abline(h=seq(from=-4,to=8,by=2),lty=2,col='grey')
abline(h=0)
# add tick every day
abline(v=dates_perDay,lty=2,col='lightblue')
# replot Darcy velocity on top
lines(x = data_flux_raw$dates,y = data_flux_raw$darcy)

# start plotting heat exchanges on 2nd axis
par(new=T)
plot(x = data_meas$dates,y=data_flux_raw$cond,lty=2,
     xaxt='n',yaxt='n',xlab='',ylab='',
     xlim=c(time1,time2),
     ylim=c(min(c(data_flux_raw$cond,data_flux_raw$adv),na.rm = T),max(c(data_flux_raw$cond,data_flux_raw$adv),na.rm = T)), # similar to other y-axis to have zeros aligned
     type='l',col='salmon3')
lines(x = data_meas$dates,y=data_flux_raw$adv,col='salmon3')
axis(side = 4,col='salmon3',col.axis = 'salmon3')
abline(h=0,lty=2,col='salmon3')
mtext(text = expression('vertical heat flux [W/m'^'2'*']'),
      side = 4,line = 2.5,col='salmon3',cex=1.5)

legend('topright',lty=c(1,2),col = 'salmon3',text.col='salmon3',
       legend = c('advective','conductive'),bg = 'white')

# abline(v=data_meas$dates[c(700,2500)],col='blue',lty=2)
par(op)
dev.off()


png(file = paste0('../PLOT/flux_timeseries_in2_plot.png'),width = 1000,height = 500)
op<-par(oma=c(1,1,1,1),mar = c(1, 4, 0, 3)+0.1,cex=1,cex.main=1.5, cex.lab=1.5, cex.axis=1.2)

par(mfrow = c(2, 1))
par(mar=c(2,4,1,0)+0.1)

# plot Darcy velocity on first y-axis
plot(x = data_flux_raw$dates,y = data_flux_raw$darcy,type='l',
     xlab='',ylab='',xaxt='n',xlim=c(time1,time2),ylim=c(min(data_flux_raw$darcy,na.rm = T),max(data_flux_raw$darcy,na.rm = T)))
mtext(text = 'water velocity [cm/day]',side = 2,line = 2.5,cex=1.5)
axis.POSIXct(side = 1,x = data_flux_raw$dates,at=dates_perWeek,
             format = '%b %d %Y',cex.axis=1.25)

# add background for positive thermal gradient
for(i in 1:nrow(df_isPos)){
  datesLim_i <- data_meas$dates[as.numeric(df_isPos[i,])]
  polygon(x=c(datesLim_i,rev(datesLim_i)),
          y=c(rep(-10,2),rep(10,2)),
          col='peachpuff',border = NA
  )
}

# add horizontal grid
abline(h=seq(from=-4,to=8,by=2),lty=2,col='grey')
abline(h=0,col='grey')
# add tick every day
abline(v=dates_perDay,lty=2,col='lightblue')
# replot Darcy velocity on top
lines(x = data_flux_raw$dates,y = data_flux_raw$darcy)
# replot box around plot
box()

# start plotting heat exchanges on 2nd axis
# par(new=T)
plot(x = data_meas$dates,y=data_flux_raw$cond,lty=2,
     xaxt='n',yaxt='n',xlab='',ylab='',
     xlim=c(time1,time2),
     ylim=c(min(c(data_flux_raw$cond,data_flux_raw$adv),na.rm = T),max(c(data_flux_raw$cond,data_flux_raw$adv),na.rm = T)),
     type='l',col='salmon3')
axis.POSIXct(side = 1,x = data_meas$dates,at=dates_perWeek,
             format = '%b %d %Y',cex.axis=1.25)

# add background for positive thermal gradient
for(i in 1:nrow(df_isPos)){
  datesLim_i <- data_meas$dates[as.numeric(df_isPos[i,])]
  polygon(x=c(datesLim_i,rev(datesLim_i)),
          y=c(rep(-100,2),rep(100,2)),
          col='peachpuff',border = NA
  )
}

axis(side = 2,col='salmon3',col.axis = 'salmon3')
mtext(text = expression('vertical heat flux [W/m'^'2'*']'),
      side = 2,line = 2,col='salmon3',cex=1.5)

# add horizontal grid
abline(h=seq(from=-400,to=400,by=20),lty=2,col='grey')
abline(h=0,col='grey')
# add tick every day
abline(v=dates_perDay,lty=2,col='lightblue')
# replot Darcy velocity on top
lines(x = data_meas$dates,y=data_flux_raw$cond,col='salmon3',lty=2)
lines(x = data_meas$dates,y=data_flux_raw$adv,col='salmon3')
# replot box around plot
box()

legend('topright',lty=c(1,2),col = 'salmon3',text.col='salmon3',
       legend = c('advective heat flux','conductive heat flux'),bg = 'white')
par(op)
# abline(v=data_meas$dates[c(700,2500)],col='blue',lty=2)
dev.off()

