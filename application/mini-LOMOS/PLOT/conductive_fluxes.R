#' Plots conductive fluxes timeseries

library(ggplot2)
library(reshape2)
library(akima)
#paths
path_output <- "../GINETTE_SENSI/OUTPUT/"
path_plot <- "../PLOT/"
path_obs <- "../GINETTE_SENSI/OBS/"
path_BC <- "../GINETTE_SENSI/"
sim_name=1


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
p

ggsave("conductive_fluxes_day.png", plot = p, width = 11, height = 8)
