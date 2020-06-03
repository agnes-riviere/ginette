pressure=read.csv("S_pressure_profil_t.dat", sep = "" , header = F ,
                  na.strings ="", stringsAsFactors= F) # time z pr h

data_pres_head=data.frame(time=pressure$V1,z=pressure$V2,p_h=pressure$V3/1000/9.81)
t11700=subset(data_pres_head, time ==11700, select = c(z, p_h))
t23400=subset(data_pres_head, time ==23400, select = c(z, p_h))
t46800=subset(data_pres_head, time ==35100, select = c(z, p_h))
jpeg("Warrick_sst.jpg")
plot(t11700$p_h,t11700$z,col="darkorange",type = "l", xlab = "Pressure head (m)", ylab = "Depth (m)")
points(-0.14,-0.91,col="darkorange")
points(-0.10,-0.87,col="darkorange")
points(-0.07,-0.79,col="darkorange")
lines(t23400$p_h,t23400$z,col="deepskyblue")
points(-0.14,-1.40,col="deepskyblue")
points(-0.10,-1.35,col="deepskyblue")
points(-0.07,-1.25,col="deepskyblue")
lines(t46800$p_h,t46800$z,col="azure4")
points(-0.14,-1.83,col="azure4")
points(-0.10,-1.77,col="azure4")
points(-0.07,-1.67,col="azure4")
legend("topleft", legend=c("Sim 11 700 s", "Sim 23 400 s", "Sim 46 800 s","Warick"),
       col=c("darkorange", "deepskyblue","azure4","black"), lty = c(1,1,1,NA), pch = c(NA,NA,NA, 1), cex=0.8)

dev.off()
