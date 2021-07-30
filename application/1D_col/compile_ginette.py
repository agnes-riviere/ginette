#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
# please compile ginette in the folder 1D_col
path = os.getcwd()

print(path)
########### Setup
# /Users/mbp/Documents/my-project/python-snippets/notebook
f_param_bck=open("/home/ariviere/Programmes/ginette/application/1D_col/E_parametre_backup.dat", "r")
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
# intrinsic permeability
val_k=10e-13
# porosity
val_n=0.3
# solid thermal conductivity
# the porous media thermal conductivity is calculated bu the Woodside relationship
#l_w = 0,598	   kg.m/s-3/C
val_l=2
# Heat capacity is calculated  by the following relationship
#  c_pm= c_w r_w n + c_s r (1-n)
# density
# c_s solid specific heat capacity
# c_s=2650D+00 m2/s2/C this value is imposed to be constant but there are no way to calibrate the both parameter rho and c in the same time
#c_w=4185D+00	       m2/s2/C
#r_w=1000  kg/m3
# solid density r=val_r 
val_r=2600e+00
# Van Genuchten parameter
val_a=6.64 #m-1
val_nVG= 2.03
swres=0.01

# user-defined  parameters zone 2
# altitude of the limit between zone 1 and zone 2
thk2=-0.2
# intrinsic permeability
val_k2=10e-15
# porosity
val_n2=0.01
# Solid thermal conductivity
val_l2=4
#solid density
val_r2=2100
# Van Genuchten parameter
val_a2=6.64 #m-1
val_nVG2= 2.03
swres2=0.010


# user-defined  parameters zone 2
# altitude of the limit between zone 1 and zone 2
thk3=-0.2
# intrinsic permeability
val_k3=10e-15
# porosity
val_n3=0.01
# Solid thermal conductivity
val_l3=4
#solid density
val_r3=2100
# Van Genuchten parameter
val_a3=6.64 #m-1
val_nVG3= 2.03
swres3=0.010



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



#os.system(r"/home/ariviere/Programmes/ginette/application/1D_col/ginette")


#znew = np.arange(z_bottom+dz/2,0, dz)
#z=pd.DataFrame(np.array([Obs1, Obs2, Obs3]),
#                   columns=['z', 't_obs'])

#f = interpolate.interp1d(x, y)


# Conductivité thermique equivalent
# Saturation
# Hydraulic conductivité (m.s-1)
# capacité calorifique
# Porosité



