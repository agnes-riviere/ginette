setwd("~/Programmes/ginette/application/ZNS_rain/Clay")

col=c(brewer.pal(n = 8, name = 'Dark2'),"darkorange","deepskyblue","azure4")
pressure=read.csv("S_pressure_profil_t.dat", sep = "" , header = F ,
                  na.strings ="", stringsAsFactors= F) # time z pr h

data_pres_head=data.frame(time=pressure$V1,z=pressure$V2,p_h=pressure$V3/1000/9.81)
t900=subset(data_pres_head, time ==900, select = c(z, p_h))
t11700=subset(data_pres_head, time ==11700, select = c(z, p_h))
t23400=subset(data_pres_head, time ==23400, select = c(z, p_h))
t36900=subset(data_pres_head,time==36900,select = c(z,p_h))
t46800=subset(data_pres_head, time ==46800, select = c(z, p_h))
t58500=subset(data_pres_head, time ==58500, select = c(z, p_h))
t62100=subset(data_pres_head, time ==62100, select = c(z, p_h))
t172800=subset(data_pres_head, time ==172800, select = c(z, p_h))

plot(t900$p_h,t900$z,col=col[1],type = "l", xlab = "Pressure head (m)", ylab = "Depth (m)")
lines(t11700$p_h,t11700$z,col=col[2])
lines(t23400$p_h,t23400$z,col=col[3])
lines(t36900$p_h,t36900$z,col=col[4])
lines(t46800$p_h,t46800$z,col=col[5])
lines(t46800$p_h,t46800$z,col=col[6])
lines(t62100$p_h,t62100$z,col=col[7])

lines(t172800$p_h,t172800$z,col=col[8])
legend("topright", legend=c("Sim 900 s","Sim 11 700 s", "Sim 23 400 s", "Sim 36 900 s", "Sim 46 800 s", "Sim 58 500 s","sim 62 100 s", "sim 172 800 s"),
       col=c(col), lty = rep(1,8), pch = rep(NA,8), cex=0.8)


