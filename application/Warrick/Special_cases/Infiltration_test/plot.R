pressure=read.csv("S_pressure_profil_t.dat", sep = "" , header = F ,
                  na.strings ="", stringsAsFactors= F) # time z pr h
pressure_cst=read.csv("S_pressure_profil_t_Q_imp.dat", sep = "" , header = F ,
                  na.strings ="", stringsAsFactors= F) # time z pr h
data_pres_head=data.frame(time=pressure$V1,z=pressure$V2,p_h=pressure$V3/1000/9.81)
data_pres_head_cst=data.frame(time=pressure_cst$V1,z=pressure_cst$V2,p_h=pressure_cst$V3/1000/9.81)
t11700=subset(data_pres_head, time ==11700, select = c(z, p_h))
t23400=subset(data_pres_head, time ==23400, select = c(z, p_h))
t46800=subset(data_pres_head, time ==46800, select = c(z, p_h))
t11700_cst=subset(data_pres_head_cst, time ==11700, select = c(z, p_h))
t23400_cst=subset(data_pres_head_cst, time ==23400, select = c(z, p_h))
t46800_cst=subset(data_pres_head_cst, time ==46800, select = c(z, p_h))
jpeg("Warrick_dT.jpg")
plot(t11700$p_h,t11700$z,col="darkorange",type = "l", xlab = "Pressure head (m)", ylab = "Depth (m)")
lines(t11700_cst$p_h,t11700_cst$z,col="orange")
points(-3.4,-0.255,col="darkorange")
points(-1.7,-0.24,col="darkorange")
points(-0.8,-0.215,col="darkorange")
lines(t23400$p_h,t23400$z,col="deepskyblue")
lines(t23400_cst$p_h,t23400_cst$z,col="blue")
points(-3.4,-0.4,col="deepskyblue")
points(-1.7,-0.38,col="deepskyblue")
points(-0.8,-0.35,col="deepskyblue")
lines(t46800$p_h,t46800$z,col="azure4")
lines(t46800_cst$p_h,t46800_cst$z,col="green")
points(-3.4,-0.6,col="azure4")
points(-1.7,-0.59,col="azure4")
points(-0.8,-0.55,col="azure4")
legend("topleft", legend=c("Sim 11 700 s", "Sim 23 400 s", "Sim 46 800 s","Warick"),
       col=c("darkorange", "deepskyblue","azure4","black"), lty = c(1,1,1,NA), pch = c(NA,NA,NA, 1), cex=0.8)

dev.off()
