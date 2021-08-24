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
os.chdir('/home/ariviere/Programmes/ginette/application/1D_col/')
# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))


if os.path.isfile('ginette'):
    print ("ginette exist")
else:
    print ("ginette not exist")
    print("you must compile ginette in the current directory")
    print(" gfortran -o ginette ../../src/ginette_V2.f")


########### Setup
# /Users/mbp/Documents/my-project/python-snippets/notebook
f_param_bck=open("E_parametre_backup.dat", "r")
setup_model=f_param_bck.read()
#time step in s
dt=900
#duration of the simulation in days
nb_day=1

#in meter
z_bottom=-0.55
dz=0.01

#apply unsaturated flow and thermal 
#unsat =1 apply
#unsat=0 cancel unsaturated zone
unsat=1


#Observation in meter
Obs1=-0.1
Obs2=-0.25
Obs3=-0.40
Obs4=-0.55

nb_zone=2
# user-defined  parameters zone 1 
# intrinsic permeability
val_k=1e-11
# porosity
val_n=0.3
# solid thermal conductivity
# the porous media thermal conductivity is calculated bu the Woodside relationship
#l_w = 0,598	   kg.m/s-3/C
val_l=3
# Heat capacity is calculated  by the following relationship
#  c_pm= c_w r_w n + c_s r (1-n)
# density
# c_s solid specific heat capacity
#val_c= c_s m2/s2/C I advice to let this value constant.
# There are no way to calibrate the both parameter rho and c in the same time.
#c_w=4185D+00	       m2/s2/C
#r_w=1000  kg/m3
# solid density r=val_r 
val_r=1600
val_c=2650
# Van Genuchten parameters
val_a=6.64 #m-1
val_nVG= 2.03
val_swres=0.01

# user-defined  parameters zone 2
# altitude of the limit between zone 1 and zone 2
thk2=-0.2
# intrinsic permeability
val_k2=1e-15
# porosity
val_n2=0.01
# Solid thermal conductivity
val_l2=4
#solid density
val_r2=2100
#solid heat capacity
val_c2=2650
# Van Genuchten parameters
val_a2=6.64 #m-1
val_nVG2= 2.03
val_swres2=0.010


# user-defined  parameters zone 2
# altitude of the limit between zone 1 and zone 2
thk3=-0.6
# intrinsic permeability
val_k3=1e-12
# porosity
val_n3=0.01
# Solid thermal conductivity
val_l3=4
#solid density
val_r3=2100
#solid heat capacitylibs_gfortran = ['gfortran']
val_c3=2650
# Van Genuchten parameters
val_a3=6.64 #m-1
val_nVG3= 2.03
val_swres3=0.010



## write the parameters
nb_cell=-z_bottom/dz
cell1=-Obs1/dz
cell2=-Obs2/dz
cell3=-Obs3/dz
cell4=-Obs4/dz


setup_model=setup_model.replace('[dt]','%06.0fD+00' % dt)
setup_model=setup_model.replace('[nb_day]','%06.0f' % nb_day)
setup_model=setup_model.replace('[z_bottom]','%6.2e' % z_bottom)
setup_model=setup_model.replace('[az]','%7.3e' % -z_bottom)
setup_model=setup_model.replace('[dz]','%6.2e' % dz)
setup_model=setup_model.replace('[nb_cell]','%05.0f' % nb_cell)
setup_model=setup_model.replace('[unsat]','%1i' % unsat)
setup_model=setup_model.replace('[cell1]','%05d' % cell1)
setup_model=setup_model.replace('[cell2]','%05d' % cell2)
setup_model=setup_model.replace('[cell3]','%05d' % cell3)
setup_model=setup_model.replace('[cell4]','%05d' % cell4)

########### Zone
f_coor=open("E_coordonnee.dat", "r")
f_zone=open("E_zone.dat", 'w')
f_paramZ_bck=open("E_zone_parameter_backup.dat", "r")
f_paramZ_new = open("E_zone_parameter.dat", 'w')
f_param_new = open("E_parametre.dat", 'w')
param_zone=f_paramZ_bck.read()
coord=pd.DataFrame()
coord = pd.read_csv(f_coor, names=["id", "x", "z"], header=None, delim_whitespace=True)



param_zone=param_zone.replace('[k1]','%8.2e' % val_k)
param_zone=param_zone.replace('[n1]','%6.2f' % val_n)
param_zone=param_zone.replace('[l1]','%6.2f' % val_l)
param_zone=param_zone.replace('[c1]','%6.2f' % val_c)
param_zone=param_zone.replace('[r1]','%6.2f' % val_r)
param_zone=param_zone.replace('[a1]','%6.2f' % val_a)
param_zone=param_zone.replace('[nVG1]','%6.2f' % val_nVG)
param_zone=param_zone.replace('[swres1]','%6.2f' % val_swres)




param_zone=param_zone.replace('[k2]','%8.2e' % val_k2)
param_zone=param_zone.replace('[n2]','%6.2f' % val_n2)
param_zone=param_zone.replace('[l2]','%6.2f' % val_l2)
param_zone=param_zone.replace('[c2]','%6.2f' % val_c2)
param_zone=param_zone.replace('[r2]','%6.2f' % val_r2)
param_zone=param_zone.replace('[a2]','%6.2f' % val_a2)
param_zone=param_zone.replace('[nVG2]','%6.2f' % val_nVG2)
param_zone=param_zone.replace('[swres2]','%6.2f' % val_swres2)

param_zone=param_zone.replace('[k3]','%8.2e' % val_k3)
param_zone=param_zone.replace('[n3]','%6.2f' % val_n3)
param_zone=param_zone.replace('[l3]','%6.2f' % val_l3)
param_zone=param_zone.replace('[c3]','%6.2f' % val_c3)
param_zone=param_zone.replace('[r3]','%6.2f' % val_r3)
param_zone=param_zone.replace('[a3]','%6.2f' % val_a3)
param_zone=param_zone.replace('[nVG3]','%6.2f' % val_nVG3)
param_zone=param_zone.replace('[swres3]','%6.2f' % val_swres3)





coord['zone'] =1

if nb_zone >= 2:
    coord['zone'] = np.where(coord['z'] <= thk2, 2,coord['zone'])
    coord['zone'] = np.where(coord['z'] <= thk3, 3,coord['zone'])

display(coord)

f_param_new.write(setup_model)
f_paramZ_new.write(param_zone)
coord.zone.to_csv(f_zone, index = False, header=False)
f_zone.close()
f_paramZ_new.close()
f_paramZ_bck.close()
f_param_bck.close()
f_param_new.close()
f_coor.close()

subprocess.call(["./ginette"])   

