# namePoint = 'syntheticStudy'
# wd=paste0('~/workspace/Programmes/LIBS/libts/trunk/test/HOBBO_HZ1D/GINETTE_SENSI/CONVERGENCE/')
# setwd(wd)

# prjPath = paste0('../../',namePoint,'/')
# # OGSprjName = namePoint
# 
# # load data
# Sys.setenv(TZ='GMT')
# 
# Deltaz = 0.01 # in m
# Deltat = 15*60 # in s
# 
# depthVect = c(0.1,0.2,0.3,0.4)
# permeability = 1e-13

# GinetteWriteInputFiles <- function(){
  
  period = 24*60*60 # in seconds
  z=seq(from=0,to=-0.4,by=-Deltaz)
  zmid = seq(from = -Deltaz/2,
             to = z[length(z)] + Deltaz/2,
             by=-Deltaz) # midpoints of cells
  
  t_seconds = seq(from=0,to=3*period,by=Deltat)
  
  # model parameters
  # notation Luce 2013
  dH = -0.01 # delta hydraulic head -1cm (infiltrating water)
  A = 1 # top teperature amplitude, in C
  T_mu = 20 #mean temperature in celcius
  # parameter values from Edmee AmPP
  lambda_w = 0.598 # W/m/K # OK GINETTE
  lambda_s = 2.3 # W/m/K # OK GINETTE
  rho_w = 1e3 # kg / m^3 # GINETTE OK
  rho_s = 2.9e3 # kg / m^3 # DENSITE MILIEU - IDEM ??
  c_w = 4.185e3 # OK GINETTE
  c_s = 1000 # J/kg/K # OK GINETTE
  mu_w = 1e-3 # viscosite dynamique OK GINETTE
  n = 0.15 # no unit porosity # OK GINETTE
  g = 9.81 # gravity acceleration
  #
  rho_mTimesc_m = n * rho_w * c_w + (1-n) * rho_s * c_s
  lambda_m = n * lambda_w + (1-n) * lambda_s
  kappa_e = lambda_m / rho_mTimesc_m
  q = permeability * rho_w * g / mu_w * dH / z[length(z)] # Darcy
  v_t = q * rho_w * c_w / rho_mTimesc_m
  omega = 2 * pi / period
  alpha = sqrt ( v_t^4 + (4 * omega * kappa_e)^2 )
  a = 1/(2*kappa_e) * (sqrt((alpha + v_t^2)/2) - v_t)
  b = 1/(2*kappa_e) * sqrt((alpha - v_t^2)/2)
  
  #### changement de la geometrie ####
  pathParam = '../E_parametre.dat'
  pathParamTemp = '../E_parametre2temp.dat'
  file.create(pathParamTemp)
  lines = readLines(pathParam)
  for(k in 1:length(lines)){
    if(grepl('nmi=',lines[k])){
      print('remplacement du nombre de mailles actives')
      idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
      substr(lines[k], start=idxEqual+1, stop=idxEqual+5) = sprintf('%.5d',abs(z[length(z)]/Deltaz))
    }
    if(grepl('nri=',lines[k])){
      print('remplacement du nombre de lignes')
      idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
      substr(lines[k], start=idxEqual+1, stop=idxEqual+5) = sprintf('%.5d',abs(z[length(z)]/Deltaz))
    }
    if(grepl('dz=',lines[k])){
      print (k)
      print('remplacement du pas despace en z')
      idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
      dzStr = gsub('e','D',sprintf('%.2e',Deltaz))
      substr(dzStr,1,4) = sprintf('%04i',as.numeric(substr(dzStr,1,4)))
      substr(lines[k], start=idxEqual+1, stop=idxEqual+8) = dzStr
    }
    if(grepl('nmaille1=',lines[k])){ #remplacement des numeros de mailles
      print('remplacement de maille profondeur 1')
      idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
      nmaille = as.integer(depthVect[1]/Deltaz)+1
      substr(lines[k], start=idxEqual+1, stop=idxEqual+5) = sprintf('%05i',nmaille)
    }
    if(grepl('nmaille2=',lines[k])){
      print('remplacement de maille profondeur 2')
      idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
      nmaille = as.integer(depthVect[2]/Deltaz)+1
      substr(lines[k], start=idxEqual+1, stop=idxEqual+5) = sprintf('%05i',nmaille)
    }
    if(grepl('nmaille3=',lines[k])){
      print('remplacement de maille profondeur 3')
      idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
      nmaille = as.integer(depthVect[3]/Deltaz)+1
      substr(lines[k], start=idxEqual+1, stop=idxEqual+5) = sprintf('%05i',nmaille)
    }
    write.table(lines[k],file=pathParamTemp,append=T,quote=F,row.names=F,col.names=F)
  }
  file.copy(from=pathParamTemp,to=pathParam,overwrite = T)
  file.remove(pathParamTemp)
  
  #### changement des conditions initiales ####

  # initial conditions - run steady state
  # change Ginette's parameters for steady state
  system('./prepGinette_ss.sh')
  # setup boundary conditions for Ginette steady state
  presBottom_ss = (dH - z[length(z)]) * rho_w * g
  tempBC_ss = T_mu + A * exp(-a*(-z[c(1,length(z))])) * cos(- b * (-z[c(1,length(z))]))
  pathCL = '../E_cdt_aux_limites.dat'
  pathCLTemp = '../E_cdt_aux_limitestemp.dat'
  file.create(pathCLTemp)
  lines = readLines(pathCL)
  for(k in 1:length(lines)){
#     if(grepl('valcl_bas=',lines[k])){ # IL FAUT METTRE EN CHARGE!! cf dessous
#       print('replace pressure condition bottom')
#       print(k)
#       idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
#       presBottom_ssIntStr = as.character(as.integer(presBottom_ss))
#       if(nchar(presBottom_ssIntStr)<=5){
#         presBottom_ssStr = paste0(paste(rep('0',5-nchar(presBottom_ssIntStr)),collapse = ''),
#                                      presBottom_ssIntStr,
#                                      'd+00')
#       }else{stop('error in defining number format in pressure condition bottom')}
#       substr(lines[k], start=idxEqual+1, stop=idxEqual+9) = presBottom_ssStr
#     }
    if(grepl('valcl_bas=',lines[k])){ # meme chose en charge
      print('replace pressure condition bottom')
      print(k)
      idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
      chargeBottom_ssIntStr = as.character(as.integer(dH*1000))
      if(nchar(chargeBottom_ssIntStr)<=5){
        if(substr(chargeBottom_ssIntStr[1],1,1) %in% c('+','-')){
          chargeBottom_ssStr = paste0(substr(chargeBottom_ssIntStr[1],1,1),
                                      paste(rep('0',5-nchar(chargeBottom_ssIntStr)),collapse = ''),
                                      substr(chargeBottom_ssIntStr,2,nchar(chargeBottom_ssIntStr)),
                                      'd-03')
        }else{
          chargeBottom_ssStr = paste0(paste(rep('0',5-nchar(chargeBottom_ssIntStr)),collapse = ''),
                                      chargeBottom_ssIntStr,
                                      'd-03')
        }
      }else{stop('error in defining number format in pressure condition bottom')}
      substr(lines[k], start=idxEqual+1, stop=idxEqual+9) = chargeBottom_ssStr
    }
    if(grepl('valclt_haut=',lines[k])){
      print('temperature condition top')
      print(k)
      idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
      tempTop_ssIntStr = as.character(as.integer(tempBC_ss[1]*1000))
      if(nchar(tempTop_ssIntStr)<=5){
        tempTop_ssStr = paste0(paste(rep('0',5-nchar(tempTop_ssIntStr)),collapse = ''),
                                  tempTop_ssIntStr,
                                  'd-03')
      }else{stop('error in defining number format in temperature condition top')}
      substr(lines[k], start=idxEqual+1, stop=idxEqual+9) = tempTop_ssStr
    } #00004D-01
    if(grepl('valclt_bas=',lines[k])){
      print('temperature condition bottom')
      print(k)
      idxEqual = gregexpr('=',lines[k])[[1]][1] # find location of equal
      tempBot_ssIntStr = as.character(as.integer(tempBC_ss[2]*1000))
      if(nchar(tempBot_ssIntStr)<=5){
        tempBot_ssStr = paste0(paste(rep('0',5-nchar(tempBot_ssIntStr)),collapse = ''),
                               tempBot_ssIntStr,
                               'd-03')
      }else{stop('error in defining number format in temperature condition bottom')}
      substr(lines[k], start=idxEqual+1, stop=idxEqual+9) = tempBot_ssStr
    }
    write.table(lines[k],file=pathCLTemp,append=T,quote=F,row.names=F,col.names=F)
  }
  file.copy(from=pathCLTemp,to=pathCL,overwrite = T)
  file.remove(pathCLTemp)

  setwd('../')
  system('gfortran -o ginette ginette_V1_10.f')
  system('./ginette')
  setwd('CONVERGENCE/')
  # change Ginette's parameters for transient simulation
  system('./prepGinette_transient.sh')

  # initial conditions temperature
  tempIC =   T_mu + A * exp(-a*(-zmid)) * cos(omega * 0 - b * (-zmid))
  write.table(x = sprintf('%2.6f',tempIC),file = '../E_temperature_initiale.dat',
              row.names = F,col.names=F,quote = F)
  
  #### changement des conditions au bord ####
  
  # calculate boundary conditions
  # boundary condition pressure
  # top -> 0
  # bottom -> delta hydraulic head -1cm (infiltrating water)
  presBC = cbind(rep('0',length.out=length(t_seconds)),
                 rep(sprintf('%2.6f',dH),length.out=length(t_seconds)))
  write.table(x = presBC,file = '../E_charge_t.dat',
            row.names = F,col.names=F,quote = F)
  
  # boundary condition temperature
  tempBC = matrix('',ncol=2,nrow=length(t_seconds))
  tempBC[,1] = sprintf('%2.6f',T_mu + A * cos(omega * t_seconds))
  tempBC[,2] = sprintf('%2.6f',T_mu + A * exp(-a*(-z[length(z)])) * 
    cos(omega * t_seconds - b * (-z[length(z)])))
  write.table(x = tempBC,file = '../E_temp_t.dat',sep='\t',
            row.names = F,col.names=F,quote = F)

# }
