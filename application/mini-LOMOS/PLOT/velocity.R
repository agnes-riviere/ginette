#' Plots velocity timeseries

library(ggplot2)
library(reshape2)
library(akima)
path_output <- "../GINETTE_SENSI/OUTPUT/"
path_plot <- "./"

sim_name=1


setwd(path_plot)
D_sim<-read.table(paste0(path_output,"Sim_velocity_profil_t_",sim_name,".dat"),header = FALSE)

#interpolation
 x=D_sim[,1]/86400
y=D_sim[,2]
z=D_sim[,3]
fld <- with(D_sim, interp(x = x, y = y, z = z))


# prepare data in long format
D_sim <- melt(fld$z, na.rm = TRUE)
names(D_sim) <- c("x", "y", "velocity")
D_sim$time <- fld$x[D_sim$x]
D_sim$depth <- fld$y[D_sim$y]

p <- ggplot(data = D_sim, aes(x = time, y = depth, z = velocity)) +
  geom_tile(aes(fill = velocity)) +
  stat_contour(colour = "black") +
  xlab("Time (days)") +
  ylab("Depth (m)") + 
   scale_fill_gradientn(colours = terrain.colors(10),name="Darcy velocity") +                              
  theme_bw()+
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 10)) 
p
ggsave("velocity_profile.png", plot = p, width = 11, height = 8)