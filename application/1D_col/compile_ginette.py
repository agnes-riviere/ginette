#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt

path = os.getcwd()

print(path)
########### Setup
# /Users/mbp/Documents/my-project/python-snippets/notebook
f_param_bck=open("/home/ariviere/Programmes/ginette/application/1D_col/E_parametre_backup.dat", "r")
f_param_new = open("/home/ariviere/Programmes/ginette/application/1D_col/E_parametre.dat", 'w')
setup_model=f_param_bck.read()
#time step in s
dt=900
#duration of the simulation in days
nb_day=1

#in meter
z_bottom=-0.27
dz=0.01

#Observation
Obs1=[-0.1, 15]
Obs1=[-0.1, 15]
Obs2=[-0.15, 12]
Obs3=[-0.2, 11]

nb_zone=2
# user-defined  parameters zone 1 
val_k=10e-6
val_n=0.3
val_l=2
val_r=2600e+00
alpha=
n=
swres

# c_s capacite calorifique specifique du solide
# c_s  impose a 1000D+00 J kg-1 k-1
#  c_pm= c_w r_w n + c_s r (1-n)  
# user-defined  parameters zone 2     
thk2=-0.2
val_k2=10e-8
val_n2=0.01
val_l2=4
val_r2=2100


nb_cell=-z_bottom/dz
   
     
setup_model=setup_model.replace('[dt]','%06.0fD+00' % dt)
setup_model=setup_model.replace('[nb_day]','%06.0f' % nb_day)
setup_model=setup_model.replace('[z_bottom]','%6.2e' % z_bottom)
setup_model=setup_model.replace('[dz]','%6.2e' % dz)
setup_model=setup_model.replace('[nb_cell]','%05.0f' % nb_cell)



########### Zone
f_coor=open("/home/ariviere/Programmes/ginette/application/1D_col/E_coordonnee.dat", "r")
f_zone=open("/home/ariviere/Programmes/ginette/application/1D_col/E_zone.dat", 'w')
f_paramZ_bck=open("/home/ariviere/Programmes/ginette/application/1D_col/E_zone_parameter_backup.dat", "r")
f_paramZ_new = open("/home/ariviere/Programmes/ginette/application/1D_col/E_zone_parameter.dat", 'w')
param_zone=f_paramZ_bck.read()
coord = pd.read_csv(f_coor, names=["id", "x", "z"], header=None, delim_whitespace=True)


  
param_zone=param_zone.replace('[k1]','%8.2e' % val_k)
param_zone=param_zone.replace('[n1]','%6.2f' % val_n)
param_zone=param_zone.replace('[l1]','%6.2f' % val_l)
param_zone=param_zone.replace('[r1]','%6.2f' % val_r)




param_zone=param_zone.replace('[k2]','%8.2e' % val_k2)
param_zone=param_zone.replace('[n2]','%6.2f' % val_n2)
param_zone=param_zone.replace('[l2]','%6.2f' % val_l2)
param_zone=param_zone.replace('[r2]','%6.2f' % val_r2)



  
     
if nb_zone == 1 :
       coord['zone'] =1   

elif nb_zone >= 2:
    coord['zone'] = [1 if s >=thk2 else 2 for s in coord['z']]   





f_param_new.write(setup_model)
f_paramZ_new.write(param_zone)
coord.zone.to_csv(f_zone, index = False)
f_paramZ_new.close()
f_paramZ_bck.close()



os.system(r"/home/ariviere/Programmes/ginette/application/1D_col/ginette")


znew = np.arange(z_bottom+dz/2,0, dz)
z=pd.DataFrame(np.array([Obs1, Obs2, Obs3]),
                   columns=['z', 't_obs'])

f = interpolate.interp1d(x, y)


# Conductivité thermique equivalent
# Saturation
# Hydraulic conductivité (m.s-1)
# capacité calorifique
# Porosité



