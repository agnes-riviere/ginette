library(RColorBrewer)
setwd("~/Programmes/ginette/application/ZNS_rain/sand")

col=c(brewer.pal(n = 8, name = 'Dark2'),"darkorange","deepskyblue","azure4")

saturation=read.csv("S_saturation_profil_t.dat", sep = "" , header = F ,
                    na.strings ="", stringsAsFactors= F) # time z pr h
xlim=c(0.1,1)
data_sat_head=data.frame(time=saturation$V1,z=saturation$V2,sat=saturation$V3)
t900=subset(data_sat_head, time ==900, select = c(z, sat))
t11700=subset(data_sat_head, time ==11700, select = c(z, sat))
t23400=subset(data_sat_head, time ==23400, select = c(z, sat))
t36900=subset(data_sat_head,time==36900,select = c(z,sat))
t46800=subset(data_sat_head, time ==46800, select = c(z, sat))
t55800=subset(data_sat_head, time ==55800, select = c(z, sat))
t58500=subset(data_sat_head, time ==58500, select = c(z, sat))
t62100=subset(data_sat_head, time ==62100, select = c(z, sat))
t86400=subset(data_sat_head, time ==86400, select = c(z, sat))
t172800=subset(data_sat_head, time ==172800, select = c(z, sat))



plot(t900$sat,t900$z,col=col[1],type = "l", xlab = "saturation ", ylab = "Depth (m)",xlim=xlim)
lines(t11700$sat,t11700$z,col=col[2])
lines(t23400$sat,t23400$z,col=col[3])
lines(t36900$sat,t36900$z,col=col[4])
lines(t46800$sat,t46800$z,col=col[5])
lines(t55800$sat,t55800$z,col=col[6])
lines(t58500$sat,t58500$z,col=col[7])
lines(t62100$sat,t62100$z,col=col[8])
lines(t86400$sat,t86400$z,col=col[9])
lines(t172800$sat,t172800$z,col=col[10])
days <- seq(1, 40, by = 1)
for (i in days) {
  t_time <- subset(data_sat_head, time == 86400 /2* i, select = c(z, sat))
  
  # Plot the subsetted data
  lines(t_time$sat, t_time$z,col='black')
  
}

legend("bottomleft", legend=c("Sim 900 s","Sim 11 700 s", "Sim 23 400 s", "Sim 46 800 s"),
       col=c("darkorange", 'red',"deepskyblue","azure4"), lty = c(1,1,1,1), pch = c(NA,NA,NA,NA), cex=0.8)



