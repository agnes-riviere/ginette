#!/usr/bin/env python

import os
import numpy as np
from pathlib import Path
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
from IPython.display import display
import subprocess
libs_gfortran = ['gfortran']
# please compile ginette in the folder 1D_col
# path of the 1D_col directory
os.chdir('/home/ariviere/Programmes/ginette/application/ZNS-1D/')
# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))


if os.path.isfile('ginette'):
    print ("ginette exist")
else:
    print ("ginette not exist")
    print("you must compile ginette in the current directory")
    print(" gfortran -o ginette ../../src/ginette_V2.f")
    subprocess.call(["gfortran","-o","ginette","../../src/ginette_V2.f"])  #creat


########### Setup

#time step in s
dt=900
#duration of the simulation in days
nb_day=100

# state
## 0 steady state
# 1 transient state (dynamic state)
state=0

#in meter
z_top=40
dz=0.033

#apply unsaturated flow and thermal 
#unsat =1 apply
#unsat=0 cancel unsaturated zone
unsat=1


nb_zone=1
# user-defined  parameters zone 1 
# intrinsic permeability [m2]  k=K*mu/(rho*g)
# K hydraulic conductivity [m.s-1]
# mu viscosity [Pa.s]
# rho density [kg.m-3]
# g gravity  9.81 [m2.s-1]
val_k=3.33333333333333e-15
# porosity
val_n=0.38 # \Phi
# solid grain density rho_s=val_r  [kg.m-3]
val_r=2578
# Van Genuchten parameters
val_a=2.00000 #m-1 alpha_vg
val_nVG= 1.23  # n_vg
val_swres=0.26 # S_wr



#Boundary conditions hydraulic head h=P/rho g + Z
top=15 
bot=15


## write the parameters
f_param_bck=open("E_parametre_backup.dat", "r")
setup_model=f_param_bck.read()

nb_cell=z_top/dz


setup_model=setup_model.replace('[dt]','%06.0fD+00' % dt)
setup_model=setup_model.replace('[state]','%1i' % state)
setup_model=setup_model.replace('[nb_day]','%06.0f' % nb_day)
setup_model=setup_model.replace('[z_top]', '%7.3e' % z_top)
setup_model=setup_model.replace('[az]','%7.3e' % z_top)
setup_model=setup_model.replace('[dz]','%6.2e' % dz)
setup_model=setup_model.replace('[nb_cell]','%05.0f' % nb_cell)
setup_model=setup_model.replace('[unsat]','%1i' % unsat)

## write the boundary conditions
f_bc_bck=open("E_cdt_aux_limites_backup.dat", "r")
bc_model=f_bc_bck.read()

bc_model=bc_model.replace('[top]', '%08.0fD+00' % top)
bc_model=bc_model.replace('[bot]','%08.0fD+00' % bot)

########### Zone
f_coor=open("E_coordonnee.dat", "r")
f_zone=open("E_zone.dat", 'w')
f_paramZ_bck=open("E_zone_parameter_backup.dat", "r")
f_paramZ_new = open("E_zone_parameter.dat", 'w')
f_param_new = open("E_parametre.dat", 'w')
f_bc_new = open("E_cdt_aux_limites.dat", 'w')
param_zone=f_paramZ_bck.read()
coord=pd.DataFrame()    
coord = pd.read_csv(f_coor, names=["id", "x", "z"], header=None, delim_whitespace=True)



param_zone=param_zone.replace('[k1]','%8.2e' % val_k)
param_zone=param_zone.replace('[n1]','%6.2f' % val_n)
param_zone=param_zone.replace('[r1]','%6.2f' % val_r)
param_zone=param_zone.replace('[a1]','%8.2e' % val_a)
param_zone=param_zone.replace('[nVG1]','%6.2f' % val_nVG)
param_zone=param_zone.replace('[swres1]','%6.2f' % val_swres)









coord['zone'] =1


display(coord)

f_param_new.write(setup_model)
f_paramZ_new.write(param_zone)
f_bc_new.write(bc_model)
coord.zone.to_csv(f_zone, index = False, header=False)
f_zone.close()
f_paramZ_new.close()
f_paramZ_bck.close()
f_param_bck.close()
f_param_new.close()
f_bc_bck.close()
f_bc_new.close()
f_coor.close()

subprocess.call(["./ginette"])   

#saturation_profile = pd.read_csv('S_saturation_profil_t.dat',sep=' ', index_col=0)


#plt.figure()
#plt.style.use('seaborn')
#Tair = df['Orly']
#Teau = df['Temp']
#plt.scatter(Tair, Teau, s=10, alpha=1, color='mediumblue',marker='.')
#plt.xlabel('Tair (°c)')
#plt.ylabel('Teau (°c)')
#df.plot(x='Orly', y='Temp', figsize=(19,4))
#plt.title('Liens Tair et Teau pour la Seine à Paris')