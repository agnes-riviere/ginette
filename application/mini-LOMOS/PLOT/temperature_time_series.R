library(ggplot2)
library(stringr)
library(RColorBrewer)
#paths
path_output <- "../GINETTE_SENSI/OUTPUT/"
path_plot <- "./"
path_obs <- "../GINETTE_SENSI/OBS/"

setwd(path_plot)
sim_name = 1
date_bg = "03/05/2019 00:00:00"
date_bg = as.POSIXct(date_bg, '%d/%m/%Y %H:%M', tz = 'GMT')

files_obs <- list.files(path = path_obs, pattern = 'Obs')
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
Melted_obs_t$time = Melted_obs_t$time + date_bg



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
Melted_sim_t$time = Melted_sim_t$time + date_bg
#Graph
## Define color
colpal <-  brewer.pal(3, "Dark2")

##plot
g_temp_ts <-
  ggplot() +
  geom_line(data = Melted_sim_t,
            mapping = aes(x = time, y = Temperature, color = Depth))  +
  geom_line(
    data = Melted_obs_t,
    mapping = aes(x = time, y = Temperature, color = Depth),
    linetype = "dashed"
  ) +
  scale_color_manual(values = c(colpal, colpal)) +
  labs(x = "", y =  expression('Temperature ('*~degree*C*')'), color = "depth") +
  scale_x_datetime(date_labels = " %d %b") +
  theme_bw() 
g_temp_ts
