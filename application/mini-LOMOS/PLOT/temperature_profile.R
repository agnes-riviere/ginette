#' Plots temperature profile timeseries

library(ggplot2)
library(reshape2)
library(akima)

path_output <- '/home/ariviere/Programmes/ginette/application/mini-LOMOS/GINETTE_SENSI/OUTPUT/'
 path_plot <- "./"

sim_name=1

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
p

ggsave("./temperature_profile.png", plot = p, width = 11, height = 8)


#library(directlabels)
#direct.label(p, list("far.from.others.borders", "calc.boxes", "enlarge.box", 
#                     hjust = 1, vjust = 1, box.color = NA, fill = "transparent", "draw.rects"))
