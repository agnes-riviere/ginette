library(ggplot2)
library(stringr)
library(RColorBrewer)
library(data.table)
library(reshape2)
library(akima)
library(cowplot)
#paths
path_output <- "../GINETTE_SENSI/OUTPUT/"
path_plot <- "../PLOT/"
path_obs <- "../GINETTE_SENSI/OBS/"
path_BC <- "../GINETTE_SENSI/"


temperature_ts <-function(sim_name = 1,date_begin= "11/04/2022 00:00:00") {
  date_begin = as.POSIXct(date_begin, '%d/%m/%Y %H:%M', tz = 'GMT')
  
  files_obs <- list.files(path = path_obs, pattern = 'Obs')
  files_BC <- list.files(path = path_BC, pattern = 'E_temp_t.dat')
  files_output <-
    list.files(path = path_output, pattern = 'temperature_maille')
  pos_end = str_locate(files_output, ".dat")
  pos_bg = str_locate(files_output, "Sim_temperature_maille")
  
  nPT100 = length(files_obs)
  # Depth PT100
  Inversion_PT100 = "inversion_PT100.COMM"
  depth_PT100 = read.csv(paste0("../", Inversion_PT100),
                         sep = " ",
                         header = FALSE) #écart entre les PT100 en cm
  depth_PT100[1, 1] = -depth_PT100[1, 1]
  if (nPT100 >= 2) {
    for (i in seq(2, nPT100 + 1)) {
      depth_PT100[i, 1] = depth_PT100[i - 1, 1] - depth_PT100[i, 1]
    }
    
  }
  depth_PT100 = depth_PT100 / 100
  
  
  #obs dataframe
  df_obs_t = read.table(paste0(path_obs, files_obs[1]),
                        header = FALSE)
  
  
  
  if (nPT100 >= 2) {
    for (i in seq(2, nPT100)) {
      a <-
        read.table(paste0(path_obs, "Obs_temperature_maille", i, "_t.dat"))
      df_obs_t <- cbind(df_obs_t, a$V2)
    }
  }
  if (nPT100 < 2)
    colnames(df_obs_t) =  c("time", paste0(depth_PT100$V1[1], " m"))
  if (nPT100 == 2)
    colnames(df_obs_t) =  c("time",
                            paste0(depth_PT100$V1[1], " m"),
                            paste0(depth_PT100$V1[2], " m"))
  if (nPT100 >= 3)
    colnames(df_obs_t) = c(
      "time",
      paste0(depth_PT100$V1[1], " m"),
      paste0(depth_PT100$V1[2], " m"),
      paste0(depth_PT100$V1[3], " m")
    )
  Melted_obs_t <- reshape2::melt(df_obs_t, id.var = "time")
  colnames(Melted_obs_t) = c("time", "Depth", "Temperature")
  Melted_obs_t$time = Melted_obs_t$time + date_begin
  
  Melted_obs_t$type='Observation' 
  
  #sim dataframe
  id_sim = as.numeric(str_sub(files_output, pos_bg[, 2] + 3,  pos_end[, 1] -
                                1))
  id = which(id_sim == sim_name)
  files_sim = files_output[id]
  pos_bg = str_locate(files_sim, "Sim_temperature_maille")
  id_pt100 = as.numeric(str_sub(files_sim, pos_bg[, 2] + 1, pos_bg[, 2] +
                                  1))
  
  df_sim_t <-  read.table(paste0(path_output, files_sim[1]),
                          header = FALSE) 
  df_BC_t = fread(paste0(path_BC, files_BC[1]),
                  header = FALSE)
  df_BC_t =df_BC_t[1:dim(df_sim_t)[1],]
  
  
  
  
  df_BC_t =data.frame(time=df_sim_t$V1,Stream=df_BC_t$V1,Bottom=df_BC_t$V2)
  
  
  if (nPT100 >= 2) {
    for (i in seq(2, nPT100)) {
      a <-  read.table(paste0(path_output, files_sim[i]),
                       header = FALSE)
      df_sim_t = cbind(df_sim_t, a$V2)
    }
  }
  if (nPT100 < 2) {
    colnames(df_sim_t) =  c("time", paste0(depth_PT100$V1[1], " m"))}
  if (nPT100 == 2) {
    colnames(df_sim_t) =  c("time",
                            paste0(depth_PT100$V1[1], " m"),
                            paste0(depth_PT100$V1[2], " m"))}
  if (nPT100 >= 3){
    colnames(df_sim_t) = c(
      "time",
      paste0(depth_PT100$V1[1], " m"),
      paste0(depth_PT100$V1[2], " m"),
      paste0(depth_PT100$V1[3], " m")
    )}
  
  Melted_sim_t <- reshape2::melt(df_sim_t, id.var = "time")
  colnames(Melted_sim_t) = c("time", "Depth", "Temperature")
  Melted_sim_t$time = Melted_sim_t$time + date_begin
  Melted_sim_t$type='Simulation'
  
  df_BC_t$time = df_BC_t$time + date_begin
  Melted_BC_t <- reshape2::melt(df_BC_t, id.var = "time")
  colnames(Melted_BC_t)  = c("time", "Depth", "Temperature")
  Melted_BC_t$type='Observation'
  #Graph
  ## Define color
  colpal <- c( brewer.pal(3, "Dark2"),"blue","red")
  if(nPT100<3){colpal<-colpal[-nPT100]}
  if(nPT100<2){colpal<-colpal[-nPT100]}
  #display.brewer.all(colorblindFriendly = TRUE)
  #display.brewer.pal(n = 8, name = 'Set2')
  #colpal <- rainbow(5)
  ##plot
  g_temp_ts <-
    ggplot()  +
    #  geom_line(data = Melted_BC_t,
    #            mapping = aes(x = time, y = Temperature, color = Depth,linetype = type))   +
    geom_line(data = Melted_sim_t,
              mapping = aes(x = time, y = Temperature, color = Depth,
                            linetype = type))  + 
    geom_line(
      data = Melted_obs_t,
      mapping = aes(x = time, y = Temperature, color = Depth,
                    linetype = type)  ) +
    scale_color_manual(values = c(colpal, colpal)) +
    labs(x = "", y =  expression('Temperature ('*~degree*C*')'), color = "depth",linetype="") +
    scale_x_datetime(date_labels = " %d %b") +
    theme_bw() 
#uncomment to print
  #  g_temp_ts
  
  ggsave("../PLOT/temperature_time_series.png", plot = g_temp_ts, width = 11, height = 8)
}


#' Plots temperature profile timeseries
temperature_profile <-function(sim_name = 1) {

# 
# setwd(path_plot)
D_sim<-read.table(paste0(path_output,"Sim_temperature_profil_t_",sim_name,".dat"),header = FALSE)

#interpolation
x=D_sim[,1]/86400
y=D_sim[,2]
z=D_sim[,3]
fld <- with(D_sim, interp(x = x, y = y, z = z))


# prepare data in long format
D_sim <- reshape2 ::melt(fld$z, na.rm = TRUE)
names(D_sim) <- c("x", "y", "temperature")
D_sim$time <- fld$x[D_sim$x]
D_sim$depth <- fld$y[D_sim$y]

p <- ggplot(data = D_sim, aes(x = time, y = depth, z = temperature)) +
  geom_tile(aes(fill = temperature))+
  stat_contour(colour = "grey")  +
  xlab("Time (days)") +
  ylab("Depth (m)") +
  scale_fill_distiller(palette = "Spectral",name = "Temperature (°C)") +
  theme_bw()+
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10))
#Uncomment to print
# p

ggsave("../PLOT/temperature_profile_day.png", plot = p, width = 11, height = 8)


#library(directlabels)
#direct.label(p, list("far.from.others.borders", "calc.boxes", "enlarge.box", 
#                     hjust = 1, vjust = 1, box.color = NA, fill = "transparent", "draw.rects"))

}


conductive_fluxes<-function(sim_name=1){
#' Plots conductive fluxes timeseries


setwd(path_plot)
D_sim<-read.table(paste0(path_output,"Sim_heat_flux_profil_t_",sim_name,".dat"),header = FALSE)

#interpolation
x=D_sim[,1]/86400
y=D_sim[,2]
z=D_sim[,3]
fld <- with(D_sim, interp(x = x, y = y, z = z))


# prepare data in long format
D_sim <- melt(fld$z, na.rm = TRUE)
names(D_sim) <- c("x", "y", "conduction")
D_sim$time <- fld$x[D_sim$x]
D_sim$depth <- fld$y[D_sim$y]
bl <- colorRampPalette(c("lightskyblue","royalblue","navy"))(200)                      
re <- colorRampPalette(c("darkred", "red2","mistyrose"))(200)
p <- ggplot(data = D_sim, aes(x = time, y = depth, z = conduction)) +
  geom_tile(aes(fill = conduction)) +
  stat_contour(colour = "black") +
  xlab("Time (days)") +
  ylab("Depth (m)") +
  scale_fill_gradientn(colours=c(re,"white", bl), na.value = "grey98",
                       limits = range(pretty(c(min(fld$z), max(fld$z)))),breaks=pretty(c(min(fld$z), max(fld$z))))   +                              
  theme_bw()+
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) 
#Uncomment to print
#p

ggsave("../PLOT/conductive_fluxes_day.png", plot = p, width = 11, height = 8)
}



advective_fluxes<-function(sim_name=1) {
#' Plots advective fluxes timeseries


setwd(path_plot)
D_sim<-read.table(paste0(path_output,"Sim_heat_flux_profil_t_",sim_name,".dat"),header = FALSE)

#interpolation
x=D_sim[,1]/86400
y=D_sim[,2]
z=D_sim[,4]
fld <- with(D_sim, interp(x = x, y = y, z = z))


# prepare data in long format
D_sim <- melt(fld$z, na.rm = TRUE)
names(D_sim) <- c("x", "y", "advection")
D_sim$time <- fld$x[D_sim$x]
D_sim$depth <- fld$y[D_sim$y]
bl <- colorRampPalette(c("lightskyblue","royalblue","navy"))(200)                      
re <- colorRampPalette(c("darkred", "red2","mistyrose"))(200)
p <- ggplot(data = D_sim, aes(x = time, y = depth, z = advection)) +
  geom_tile(aes(fill = advection)) +
  stat_contour(colour = "black") +
  xlab("Time (days)") +
  ylab("Depth (m)") +
  scale_fill_gradientn(colours=c(re,"white", bl), na.value = "grey98",
                       limits = range(pretty(c(min(fld$z), max(fld$z)))),breaks=pretty(c(min(fld$z), max(fld$z))))  +                              
  theme_bw()+
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) 
#Uncomment to print
#p
ggsave("../PLOT/advective_fluxes_day.png", plot = p, width = 11, height = 8)
}

flux_ts<-function(  namePointT='T3_Point3Nelly_14_04_22.csv',namePointP='P2_Point3Nelly_14_04_22.csv',sim_name=1,date_begin='14/04/2022 17:45:00')
  {
  #paths
  path_output <- "../GINETTE_SENSI/OUTPUT/"
  path_plot <- "../PLOT/"
  path_obs <- "../GINETTE_SENSI/OBS/"
  path_BC <- "../GINETTE_SENSI/"

  pathFile_velocity =paste0(path_output,'S_vitesse_nmaille2_hb_',sim_name,'.dat')
  
  date_begin=as.POSIXct(date_begin,tz = 'GMT',format='%d/%m/%Y %H:%M:%S')
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
  
  time1 = date_begin+86400
  time2 = max(data_flux_raw$dates,na.rm=T)
  # -- read files to get dates --
  
  pathMeas = paste0(path_BC,'../',namePointT)
  data_meas <- fread(file = pathMeas,header = T,sep = ',')
  data_meas=data_meas[, Filter(function(x) any(!is.na(x)), .SD)]
  colnames(data_meas)[1]<-'id'
  colnames(data_meas)[2]<-'dates'
  if( ncol(data_meas)>=3) {colnames(data_meas)[3]<-'T_depth_1_C'}
  if( ncol(data_meas)>=4)  colnames(data_meas)[4]<-'T_depth_2_C'
  if( ncol(data_meas)>=5)  colnames(data_meas)[5]<-'T_depth_3_C'
  if( ncol(data_meas)>=6)  colnames(data_meas)[6]<-'T_depth_4_C'
  data_meas$dates <- as.POSIXct(x = data_meas$dates,tz = 'GMT',format='%d/%m/%Y %H:%M:%S')
  data_meas=data_meas[data_meas$dates >= time1&data_meas$dates <= time2, ]
  pathMeasP = paste0(path_BC,'../',namePointP)
  data_measP <- fread(file = pathMeasP,header = T,sep = ',')
  data_measP=data_measP[, Filter(function(x) any(!is.na(x)), .SD)]
  colnames(data_measP)=c('id','dates','dh','T_stream_C')
  data_measP$dates <- as.POSIXct(x = data_measP$dates,tz = 'GMT',format='%d/%m/%Y %H:%M:%S')
  
  
  time1 = max(date_begin+86400,min(data_measP$dates))
  time2 = min(c(max(data_flux_raw$dates,na.rm=T),max(data_meas$dates,na.rm=T),max(data_measP$dates,na.rm=T)))
  
  date_end=min(time2,max(data_meas$dates),max(data_measP$dates))
  dates_perDay <- seq(from=time1,
                      to=date_end,
                      by="day")
  dates_perWeek <- dates_perDay[seq(from=1,to=length(dates_perDay),by=7)]
  data_measP=data_measP[data_measP$dates >= time1 &data_measP$dates <= time2 , ]
  class(data_meas) <- class(as.data.frame(data_meas))
  class(data_measP) <- class(as.data.frame(data_measP))
  data_meas<-subset(data_meas,data_meas$dates>time1)
  data_measP<-subset(data_measP,data_measP$dates>time1)
  data_meas<-subset(data_meas,data_meas$dates<time2)
  data_measP<-subset(data_measP,data_measP$dates<time2)
  
  
  # in the plot we want to add
  # whether the temperature gradient in the river
  # is positive or negative
  
  # to do so, check the temperature timeseries
  # in the river and at the bottom of the HZ
  data_measP$T_stream_C <- as.numeric(data_measP$T_stream_C)
  data_meas[,ncol(data_meas)] <- as.numeric(data_meas[,ncol(data_meas)])
  
  plot(data_measP$T_stream_C,type='l')
  lines(data_meas[,ncol(data_meas)])
  
  # isPosTherm contains true if the thermal gradient is positive
  # (ie warmer in the stream than in the hyporheic zone)
  data_meas$isPosTherm <- ((data_measP$T_stream_C - data_meas[,ncol(data_meas)]) > 0)
  idx_pos2neg <- 1 + which(diff(data_meas$isPosTherm) == -1)
  idx_neg2pos <- 1 + which(diff(data_meas$isPosTherm) == 1)
  df_isPos <- data.frame(begin=c(idx_neg2pos), end=idx_pos2neg)
  
  for(i in 1:nrow(df_isPos)){
    polygon(x=c(df_isPos[i,],rev(df_isPos[i,])),
            y=c(rep(0,2),rep(20,2)),
            col='peachpuff',border = NA
    )
  }
  lines(data_meas$T_stream_C,type='l',col='black')
  lines(data_meas[,ncol(data_meas)])
  dev.off()
  # how to plot that now....
  # 
  
  png(file = paste0('../PLOT/flux_timeseries.png'),width = 1000,height = 500)
  op<-par(oma=c(1,1,1,1),mar = c(1, 3, 0, 3)+0.1,cex=1,cex.main=1.5, cex.lab=1.5, cex.axis=1.2)
  data_flux_raw=subset(data_flux_raw,data_flux_raw$time>86400)
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
  plot(x = data_flux_raw$dates,y=data_flux_raw$cond,lty=2,
       xaxt='n',yaxt='n',xlab='',ylab='',
       xlim=c(time1,time2),
       ylim=c(min(c(data_flux_raw$cond,data_flux_raw$adv),na.rm = T),max(c(data_flux_raw$cond,data_flux_raw$adv),na.rm = T)), # similar to other y-axis to have zeros aligned
       type='l',col='salmon3')
  lines(x = data_flux_raw$dates,y=data_flux_raw$adv,col='salmon3')
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
  plot(x = data_flux_raw$dates,y=data_flux_raw$cond,lty=2,
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
            y=c(rep(min(c(data_flux_raw$cond,data_flux_raw$adv),na.rm = T),2),rep(max(c(data_flux_raw$cond,data_flux_raw$adv),na.rm = T),2)),
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
  lines(x = data_flux_raw$dates,y=data_flux_raw$cond,col='salmon3',lty=2)
  lines(x = data_flux_raw$dates,y=data_flux_raw$adv,col='salmon3')
  # replot box around plot
  box()
  
  legend('topright',lty=c(1,2),col = 'salmon3',text.col='salmon3',
         legend = c('advective heat flux','conductive heat flux'),bg = 'white')
  par(op)
  # abline(v=data_meas$dates[c(700,2500)],col='blue',lty=2)
  dev.off()
  
  
  
}



tot_fig<-function(  namePointT='T3_Point3Nelly_14_04_22.csv',namePointP='P2_Point3Nelly_14_04_22.csv',sim_name=1,date_begin='14/04/2022 17:45:00')
{
  #paths
  path_output <- "../GINETTE_SENSI/OUTPUT/"
  path_plot <- "../PLOT/"
  path_obs <- "../GINETTE_SENSI/OBS/"
  path_BC <- "../GINETTE_SENSI/"
  
  pathFile_velocity =paste0(path_output,'S_vitesse_nmaille2_hb_',sim_name,'.dat')
  
  date_begin=as.POSIXct(date_begin,tz = 'GMT',format='%d/%m/%Y %H:%M:%S')

  # -- read files with ginette outputs --
  path_output <- "../GINETTE_SENSI/OUTPUT/"
  path_plot <- "../PLOT/"
  path_obs <- "../GINETTE_SENSI/OBS/"
  path_sim <- "../GINETTE_SENSI/"
  
  # source('utils_multiplot.R')
  
  
  
  
  # read velocity profile
  dat_velocity <- fread(paste0(path_output,"Sim_velocity_profil_t_",sim_name,".dat"),header = FALSE)
  names(dat_velocity) <- c('t_s','z','q')
  
  # read temperature profile
  dat_temperature <- fread(paste0(path_output,"Sim_temperature_profil_t_",sim_name,".dat"),header = FALSE)
  names(dat_temperature) <- c('t_s','z','T')
  
  # read heat profile
  dat_heat <- fread(paste0(path_output,"Sim_heat_flux_profil_t_",sim_name,".dat"),header = FALSE)
  D_sim <- data.frame(
    t_s=dat_heat$V1[seq(from=1,to=nrow(dat_heat),by=2)],
    z=dat_heat$V2[seq(from=1,to=nrow(dat_heat),by=2)],
    conduction=dat_heat$V3[seq(from=1,to=nrow(dat_heat),by=2)],
    advection=dat_heat$V4[seq(from=2,to=nrow(dat_heat),by=2)],
    total_heat=dat_heat$V5[seq(from=2,to=nrow(dat_heat),by=2)]
  )
  
  # group all variables in same dataframe D_sim
  D_sim <- dplyr::left_join(x=D_sim,y=dat_temperature)
  D_sim <- dplyr::left_join(x=D_sim,y=dat_velocity)
  
  # ---- also add real dates ----
  D_sim$datesPosix=D_sim$t_s+date_begin
  
  D_sim <- D_sim[!is.na(D_sim$datesPosix),]
  
  # plot
  library(ggplot2)
  
  # first plot velocity
  D_sim<-subset(D_sim,D_sim$t_s>86400)
  D_sim$q_cmPerDay <- D_sim$q * 100 * (3600*24)
  limits = c(floor(min(D_sim$q_cmPerDay)), ceiling(max(D_sim$q_cmPerDay)))
  
  if(limits[2]==1) {
    limits = c(min(D_sim$q_cmPerDay), max(D_sim$q_cmPerDay))
  }
  
  
  D_sim$q_cmPerDay <- D_sim$q * 100 * (3600*24)
  
  g_velocity <-
    ggplot(data = D_sim, aes(x = datesPosix, y = z)) +
    geom_raster(aes(fill=q_cmPerDay),interpolate = TRUE) +
    scale_x_datetime(labels = scales::date_format("%b %d %Y"),
                     date_minor_breaks = "1 day") +
    xlab(NULL) +
    ylab("Depth [m]") + 
    scale_fill_gradientn(colours = rainbow(100),
                         na.value = "grey98",
                         limits =limits,
                         name="Darcy\nvelocity\n[cm/day]") +                              
    theme_bw()+
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10)) +
    geom_text(data=NULL,aes(datesPosix[1],y=z[1],label="(a)"))
  
  ggsave("./velocity_profile.png", plot = g_velocity, width = 11, height = 8)
  
  
  # plot temperature
  g_temperature <-
    ggplot(data = D_sim, aes(x = datesPosix, y = z)) +
    geom_raster(aes(fill=T),interpolate = TRUE) +
    scale_x_datetime(labels = scales::date_format("%b %d %Y"),
                     date_minor_breaks = "1 day") +
    xlab(NULL) +
    ylab("Depth [m]") + 
    scale_fill_distiller(palette = "Spectral",name = "temperature\n[°C]") +                           
    theme_bw()+
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10))
  
  ggsave("./temperature_profile.png", plot = g_temperature, width = 11, height = 8)
  
  bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
  re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
  color_lim_ad <- c(floor(min(D_sim$advection)), ceiling(max(D_sim$advection)))
  if(color_lim_ad[1]==-1) {
    color_lim_ad <- c(min(D_sim$advection), max(D_sim$advection))
  }
  g_advection <-
    ggplot(data = D_sim, aes(x = datesPosix, y = z)) +
    geom_raster(aes(fill=advection),interpolate = TRUE) +
    scale_x_datetime(labels = scales::date_format("%b %d %Y"),
                     date_minor_breaks = "1 day") +
    xlab(NULL) +
    ylab("Depth [m]") +
    scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
                         na.value = "grey98",
                         limits = color_lim_ad,
                         name=expression('advective\nheat\nflux [W/m'^'2'*']'))  +                              
    theme_bw()+
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10))
  #+geom_text(data=NULL,aes(datesPosix[1],y=z[1]))
  
  ggsave("./advective_profile.png", plot = g_advection, width = 11, height = 8)
  
  
  bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
  re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
  color_lim_con <- c(floor(min(D_sim$conduction)), ceiling(max(D_sim$conduction)))
  if(color_lim_con[1]==-1) {
    color_lim_con <- c(min(D_sim$conduction), max(D_sim$conduction))
  }
  g_conduction <- 
    ggplot(data = D_sim, aes(x = datesPosix, y = z)) +
    geom_raster(aes(fill=conduction),interpolate = TRUE) +
    scale_x_datetime(labels = scales::date_format("%b %d %Y"),
                     date_minor_breaks = "1 day") +
    # stat_contour(colour = "black") +
    xlab(NULL) +
    ylab("Depth [m]") +
    scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
                         na.value = "grey98",
                         limits = color_lim_con,
                         name=expression('conductive\nheat\nflux [W/m'^'2'*']'))  +  
    theme_bw()+
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10)) 
  #+  geom_text(data=NULL,aes(datesPosix[1],y=z[1]))
  
  ggsave("./conductive_profile.png", plot = g_advection, width = 11, height = 8)
  
  
  bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
  re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
  color_lim_tot <- c(floor(min(D_sim$total_heat)), ceiling(max(D_sim$total_heat)))
  if(color_lim_tot[1]==-1) {
    color_lim_tot<- c(min(D_sim$total_heat), max(D_sim$total_heat))
  }
  g_sum <- 
    ggplot(data = D_sim, aes(x = datesPosix, y = z)) +
    geom_raster(aes(fill=total_heat),interpolate = TRUE) +
    scale_x_datetime(labels = scales::date_format("%b %d %Y"),
                     date_minor_breaks = "1 day") +
    # stat_contour(colour = "black") +
    xlab(NULL) +
    ylab("Depth [m]") +
    scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
                         na.value = "grey98",
                         limits = color_lim_tot,
                         name=expression('total heat\nflux [W/m'^'2'*']'))  +  
    theme_bw()+
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 10))
  ggsave("./total_heat_profile.png", plot = g_advection, width = 11, height = 8)
  #+  geom_text(data=NULL,aes(datesPosix[1],y=z[1],label="(d)"))
  
  # plot all
  # png('jolies_frises.png',width = 500,height = 700)
  # multiplot(g_temperature,g_advection,g_conduction,g_sum,cols = 1)
  # dev.off()
  
  # modification des l?gendes pour figures en fran?ais
  png("./jolie_frise_EN.png",width = 400,height=600)
  plot_grid(g_temperature + 
              scale_x_datetime(labels = scales::date_format("%d/%m/%Y"),
                               date_minor_breaks = "1 day") +
              scale_fill_distiller(palette = "Spectral",name = "température\n[°C]") +
              labs(y = 'profondeur [m]'),
            g_advection +
              scale_x_datetime(labels = scales::date_format("%d/%m/%Y"),
                               date_minor_breaks = "1 day") +
              scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
                                   na.value = "grey98",
                                   limits = color_lim_ad,
                                   name=expression('Advective\nheat\nflux  [W/m'^'2'*']')) +
              labs(y = 'profondeur [m]'),
            g_conduction +
              scale_x_datetime(labels = scales::date_format("%d/%m/%Y"),
                               date_minor_breaks = "1 day") +
              scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
                                   na.value = "grey98",
                                   limits = color_lim_con,
                                   name=expression('conductive\nheat\nflux [W/m'^'2'*']')) +
              labs(y = 'profondeur [m]'),
            g_sum +
              scale_x_datetime(labels = scales::date_format("%d/%m/%Y"),
                               date_minor_breaks = "1 day") +
              scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
                                   na.value = "grey98",
                                   limits = color_lim_tot,
                                   name=expression('Total\nheat\nflux  [W/m'^'2'*']')) +
              labs(y = 'profondeur [m]'),
            ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 12)
  dev.off()
  
  png("./jolie_frise_fr.png",width = 400,height=600)
  plot_grid(g_temperature + 
              scale_x_datetime(labels = scales::date_format("%d/%m/%Y"),
                               date_minor_breaks = "1 day") +
              scale_fill_distiller(palette = "Spectral",name = "temperature\n[°C]") +
              labs(y = 'depth [m]'),
            g_advection +
              scale_x_datetime(labels = scales::date_format("%d/%m/%Y"),
                               date_minor_breaks = "1 day") +
              scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
                                   na.value = "grey98",
                                   limits = color_lim_ad,
                                   name=expression('flux de\nchaleur\nadvectif [W/m'^'2'*']')) +
              labs(y = 'depth [m]'),
            g_conduction +
              scale_x_datetime(labels = scales::date_format("%d/%m/%Y"),
                               date_minor_breaks = "1 day") +
              scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
                                   na.value = "grey98",
                                   limits = color_lim_con,
                                   name=expression('flux de\nchaleur\nconductif [W/m'^'2'*']')) +
              labs(y = 'depth [m]'),
            g_sum +
              scale_x_datetime(labels = scales::date_format("%d/%m/%Y"),
                               date_minor_breaks = "1 day") +
              scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
                                   na.value = "grey98",
                                   limits = color_lim_tot,
                                   name=expression('flux de\nchaleur\ntotal [W/m'^'2'*']')) +
              labs(y = 'depth [m]'),
            ncol = 1,
            labels = c('A', 'B', 'C', 'D'),
            label_size = 12)
  dev.off()
  
  
}
