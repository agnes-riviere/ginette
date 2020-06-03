pressure=read.csv("S_pressure_profil_t.dat", sep = "" , header = F ,
                  na.strings ="", stringsAsFactors= F) # time z pr h

val1=28944000
val2=30758400
val3=31017600
data_pres_head=data.frame(time=pressure$V1,z=pressure$V2,p_h=pressure$V3/1000/9.81)
t11700=subset(data_pres_head, time ==val1, select = c(z, p_h))
t23400=subset(data_pres_head, time ==val2, select = c(z, p_h))
t46800=subset(data_pres_head, time ==val3, select = c(z, p_h))
write.csv(data_pres_head, "data_pres_head.csv", sep = "\t")
leg1= paste0("Simulation: ", round(val1/(3600*24)), " days")
leg2= paste0("Simulation: ", round(val2/(3600*24)), " days")
leg3= paste0("Simulation: ", round(val3/(3600*24)), " days")
legs = c(leg1,leg2,leg3)

png("ORCHIDEE.png")
plot(t11700$p_h,t11700$z,col="darkorange",type = "l", xlab = "Pressure head (m)", ylab = "Depth (m)")
lines(t23400$p_h,t23400$z,col="deepskyblue")
lines(t46800$p_h,t46800$z,col="azure4")
legend("topleft", legend=legs,
       col=c("darkorange", "deepskyblue","azure4","black"), lty = c(1,1,1,NA), pch = c(NA,NA,NA, 1), cex=0.8)

dev.off()
 