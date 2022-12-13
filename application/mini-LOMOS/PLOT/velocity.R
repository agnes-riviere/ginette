library(ggplot2)
library(reshape2)
library(data.table)
library(rgl)

namePointT <- 'T3_Point3Nelly_14_04_22.csv'
namePointP <- 'P2_Point3Nelly_14_04_22.csv'
iFile=58
isim=40
date_begin=as.POSIXct('14/04/2022 17:45:00',tz = 'GMT',format='%d/%m/%Y %H:%M:%S')
path_sim='/home/ariviere/Programmes/ginette/application/mini-LOMOS/GINETTE_SENSI/'
path_plot='/home/ariviere/Programmes/ginette/application/mini-LOMOS/PLOT/'

# setwd("C:/Data/compte_linux/Documents/in_zip/")
D_sim<-fread(paste0(path_sim,"Sim_velocity_profil_t.dat"),header = FALSE)

x=D_sim[,1]/86400
y=D_sim[,2]
z=D_sim[,3]*3600*100*24
fld <- with(D_sim, akima::interp(x = x, y = y, z = z))

# prepare data in long format
D_sim <- reshape2::melt(fld$z, na.rm = TRUE)
names(D_sim) <- c("x", "y", "velocity")
D_sim$time <- fld$x[D_sim$x]
D_sim$depth <- fld$y[D_sim$y]

p <- ggplot(data = D_sim, aes(x = time, y = depth, z = velocity)) +
  geom_tile(aes(fill = velocity)) +
  # stat_contour(colour = "black") +
  xlab("Time (days)") +
  ylab("Depth (m)") + 
  scale_fill_gradientn(colours = terrain.colors(10),name="Darcy velocity") +                              
  theme_bw()+
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) 
p
