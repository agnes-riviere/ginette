library(ggplot2)
library(reshape2)
library(data.table)
library(rgl)
library(akima)

namePointT <- 'T3_Point3Nelly_14_04_22.csv'
namePointP <- 'P2_Point3Nelly_14_04_22.csv'
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

dat_velocity$q_cmPerDay <- dat_velocity$q * 100 * (3600*24)

g_velocity <-
  ggplot(data = dat_velocity, aes(x = t_s, y = z)) +
  geom_raster(aes(fill=q_cmPerDay),interpolate = TRUE) +
  xlab('day') +
  ylab("Depth [m]") + 
  scale_fill_gradientn(colours = rainbow(100),
                       na.value = "grey98",
                       limits = c(floor(min(D_sim$q_cmPerDay)), ceiling(max(D_sim$q_cmPerDay))),
                       name="Darcy\nvelocity\n[cm/day]") +                              
  theme_bw()+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 10)) 

ggsave("./velocity_profile_day.png", plot = g_velocity, width = 11, height = 8)
