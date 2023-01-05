## plot velocity time series f

library(cowplot)
library(data.table)


isim=2

date_begin=as.POSIXct('14/04/2022 17:45:00',tz = 'GMT',format='%d/%m/%Y %H:%M:%S')
# -- read files with ginette outputs --
path_output <- "../GINETTE_SENSI/OUTPUT/"
path_plot <- "../PLOT/"
path_obs <- "../GINETTE_SENSI/OBS/"
path_sim <- "../GINETTE_SENSI/"


# read velocity profile
dat_velocity <- fread(paste0(path_output,"Sim_velocity_profil_t_",isim,".dat"),header = FALSE)
names(dat_velocity) <- c('t_s','z','q')

# read temperature profile
dat_temperature <- fread(paste0(path_output,"Sim_temperature_profil_t_",isim,".dat"),header = FALSE)
names(dat_temperature) <- c('t_s','z','T')

# read heat profile
dat_heat <- fread(paste0(path_output,"Sim_heat_flux_profil_t_",isim,".dat"),header = FALSE)
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


## change later - un peu nul mais ca marche
# j'avais oublie de changer le nombre de pas de temps dans ginette
# a faire automatiquement
# enlever les NA dans datesPosix qui correspondent a des temps en dehors des mesures
D_sim <- D_sim[!is.na(D_sim$datesPosix),]

# plot
library(ggplot2)

# first plot velocity

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
                       limits = c(floor(min(D_sim$q_cmPerDay)), ceiling(max(D_sim$q_cmPerDay))),
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

## --- plot for graphical abstract --

##png("jolie_frise_EN.png",width = 400,height=200)
##ggplot(data = D_sim, aes(x = datesPosix, y = z)) +
##  geom_raster(aes(fill=total_heat),interpolate = TRUE) +
## scale_x_datetime(labels = scales::date_format("%b %d %Y"),
##              date_minor_breaks = "1 day") +
  # stat_contour(colour = "black") +
  ##  xlab("Time") +
  ##  ylab("Depth") +
  ##scale_fill_gradientn(colours=c(rev(re),"white", rev(bl)),
  ##                       na.value = "grey98",
  ##                       limits = color_lim,
##       name=expression('total heat\nflux [W/m'^'2'*']'))  +  
  ## theme_bw()+
  ## theme(legend.title = element_text(size = 10),
  ##      legend.text = element_text(size = 10),
  ##       axis.text=element_blank(),
  ##        axis.ticks=element_blank()) 
  ##dev.off()