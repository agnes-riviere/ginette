#set wd to source file location
setwd("./")

library(hydroGOF)
library(stringr)
library(plyr)

path_obs <- "GINETTE_SENSI/OBS/"
path_output <- "GINETTE_SENSI/OUTPUT/"

files_obs <- list.files(path = path_obs, pattern = 'Obs')
files_output <- list.files(path = path_output, pattern = 'temperature_maille')

nb_PT100 <- length(files_obs)
num_PT100 <- substr(files_obs, 23, 23)

nb_simu <- as.integer(str_remove(string = str_sub(files_output[length(files_output)], 25), pattern = ".dat"))
param_simu <- read.table(file = "GINETTE_SENSI/tested_values", sep = " ", header = FALSE)

#Read data

for (i in seq_len(nb_PT100))  { 
  a <- read.csv(paste0(path_obs, files_obs)[i], sep=" ", header = FALSE)
  name <- str_remove(files_obs[i], ".dat")
  assign(name, a)
}

for (i in seq_len(nb_PT100*nb_simu))  { 
  a <- read.csv(paste0(path_output, files_output)[i], sep = " ", header = FALSE)
  name <- str_remove(files_output[i], ".dat")
  assign(name, a)
}

#Creation of a data-frame with all the observed temperatures

data_obs <- data.frame("1" = Obs_temperature_maille1_t, "2" = Obs_temperature_maille2_t[2], "3" = Obs_temperature_maille3_t[2])[-1, ]

#Creation of a data-frame with all the simulated temperatures

data_output <- Sim_temperature_maille1_1[1]

for (i in seq_len(nb_simu)) {
  for (j in seq_len(nb_PT100)) {
  a <- get(paste0(paste0("Sim_temperature_maille",j,"_"),i))[2]
  data_output <- cbind(data_output, a)
  }
}

#Creation of the matrix of the stats results of comparisons obs with sim
results <- array(NA, dim = c(nb_simu, 5*nb_PT100+1))

#name of rows
rows <- as.character(c())
for (i in seq_len(nb_simu)) {
   rows[i] <- paste0("Simu", i, " ", 
                     "k=", param_simu[i, 1],
                     " n=", param_simu[i, 2],
                     " l=", param_simu[i, 3],
                     " r=", param_simu[i, 4])
}
rownames(results) <- rows

#name of columns
columns <- as.character(c())
columns[1] <- "Sum KGE"
for (i in seq_len(nb_PT100)) {
  columns[i+1] <- paste0("KGE m", num_PT100[i])
  columns[i+1+as.integer(num_PT100[length(num_PT100)])] <- paste0("RMSE m", num_PT100[i])
  columns[i+1+2*as.integer(num_PT100[length(num_PT100)])] <- paste0("MAE m", num_PT100[i])
  columns[i+1+3*as.integer(num_PT100[length(num_PT100)])] <- paste0("COR m", num_PT100[i])
  columns[i+1+4*as.integer(num_PT100[length(num_PT100)])] <- paste0("PBIAS m", num_PT100[i])
}
colnames(results) <- columns

#Fill the matrix

if (nb_simu==1) {
  a <- c(0, 3)
} else {
  a <- seq(0, (nb_simu - 1) * 3, by = 3)
}

for (j in seq_len(nb_simu)) {
  for (i in seq_len(nb_PT100)) {
    results[j, i + 1] <- KGE(data_obs[, i + 1], data_output[, i + a[j] + 1])
    results[j, i + 4] <- rmse(data_obs[, i + 1], data_output[, i + a[j] + 1])
    results[j, i + 7] <- mae(data_obs[, i + 1], data_output[, i + a[j] + 1])
    results[j, i + 10] <- cor(data_obs[, i + 1], data_output[, i + a[j] + 1])
    results[j, i + 13] <- pbias(data_obs[, i + 1], data_output[, i + a[j] + 1])
    if (nb_PT100 == 1) {
      results[j, 1] <- results[j, 2]
    } else if (nb_PT100 == 2) {
      results[j, 1] <- results[j, 2] + results[j, 3]
    } else {
      results[j, 1] <- results[j, 2] + results[j, 3] + results[j, 4]
    }
  }
}

#Sort depending on the sum of KGE
sorted_results <- results[order(results[, 1], decreasing = TRUE),]

write.table(sorted_results, file = "SENSI/Results_Stats_sim-obs", sep = ";", col.names = TRUE, row.names = TRUE)
