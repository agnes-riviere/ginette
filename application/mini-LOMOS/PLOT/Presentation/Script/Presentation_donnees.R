#LIBRARY
#---- Interpolate timeseries to match the time discretization ----
library(stats) #for spline interpolation
library(lubridate)
library(stringr)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(data.table)
library(readr)

setwd("./")

#paths
path_data <- "../Data_MiniLomos/"
path_output <- "../Output/"

file_depth <- list.files(path = path_data, pattern = "inversion")
point_names <- str_remove(string = file_depth, pattern = "_inversion_PT100.COMM")
file_temperature <- list.files(path = path_data, pattern = "tPoint")

depth_data <- list()
for (i in 1:length(file_depth)) {
  depth_data[i] <- fread(file = paste0(path_data, file_depth[i]))
}

temperature_data <- list()
for (i in 1:length(file_temperature)) {
t_data <- as.data.frame(fread(paste0(path_data, file_temperature[i])))
t_data <- t_data[, -1]
t_data$dates <- as.POSIXct(t_data$dates, '%d/%m/%Y %H:%M', tz = 'GMT')
name <- append(c("Date", "Stream"), paste0(as.character(depth_data[[i]]), " m"), after = 2)
colnames(t_data) <- name
Melted_t_data <- reshape2::melt(t_data, id.var = "Date")
colnames(Melted_t_data) = c("time", "Depth", "Temperature")
Melted_t_data$Type <- point_names[i]
temperature_data[[i]] <- Melted_t_data
}

t_data <- fread(file = paste0(path_data, file_temperature[7]))
t_data <- t_data[, -1]
t_data$dates <- as.POSIXct(t_data$dates, '%d/%m/%y %H:%M', tz = 'GMT')
name <- append(c("Date", "Stream"), paste0(as.character(depth_data[[7]]), " m"), after = 2)
colnames(t_data) <- name
Melted_t_data <- reshape2::melt(t_data, id.var = "Date")
colnames(Melted_t_data) = c("time", "Depth", "Temperature")
Melted_t_data$Type <- point_names[7]
temperature_data[[7]] <- Melted_t_data

#Figure
colpal <-  brewer.pal(8, "Dark2")
g_temp_ts <- ggplot()
for (i in 1:length(file_temperature)) {
  g_temp_ts <- g_temp_ts + 
    geom_line(data = temperature_data[[i]],
              mapping = aes(x = time, y = Temperature, color = Type, linetype = Depth))
}

g_temp_ts <- g_temp_ts + 
  scale_color_manual(values = c(colpal)) +
  scale_x_datetime(date_labels = " %d %b %Y") +
  labs(x = "Date", y = "T (°C)", color = "Depth") +
  theme_bw()

png(paste0(path_output, "Test", ".png"), width = 720, height = 480, units = "px")
g_temp_ts
dev.off()
