pressure=read.csv("S_pressure_profil_t.dat", sep = "" , header = F ,
                  na.strings ="", stringsAsFactors= F) # time z pr h

data_pres_head=data.frame(time=pressure$V1,z=pressure$V2,p_h=pressure$V3/1000/9.81)
data_pres_head$time = round(data_pres_head$time)
t11700=subset(data_pres_head, time ==11700, select = c(z, p_h))
t23400=subset(data_pres_head, time ==23400, select = c(z, p_h))
t35100=subset(data_pres_head, time ==35100, select = c(z, p_h))
jpeg("Warrick_limestone.jpg")
plot(t11700$p_h,t11700$z,col="darkorange",type = "l", xlab = "Pressure head (m)", ylab = "Depth (m)")
points(-4.18,-0.90,col="darkorange")
points(-2.69,-0.85,col="darkorange")
points(-1.61,-0.75,col="darkorange")
lines(t23400$p_h,t23400$z,col="deepskyblue")
points(-4.18,-1.34,col="deepskyblue")
points(-2.69,-1.28,col="deepskyblue")
points(-1.61,-1.15,col="deepskyblue")
lines(t35100$p_h,t35100$z,col="azure4")
points(-4.18,-1.71,col="azure4")
points(-2.69,-1.64,col="azure4")
points(-1.61,-1.49,col="azure4")
legend("topleft", legend=c("Sim 11 700 s", "Sim 23 400 s", "Sim 35 100 s","Warick"),
       col=c("darkorange", "deepskyblue","azure4","black"), lty = c(1,1,1,NA), pch = c(NA,NA,NA, 1), cex=0.8)


dev.off()

