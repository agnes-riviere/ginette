#!/usr/bin/env python

import os
import numpy as np
from pathlib import Path
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
from IPython.display import display
import subprocess
import shutil
libs_gfortran = ['gfortran']
# please compile ginette in the folder 1D_col
# path of the 1D_col directory
path_mini_lomos='/home/ariviere/Programmes/ginette/application/mini-LOMOS/'
os.chdir(path_mini_lomos)
# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))
path_one_sim='GINETTE_SENSI'

if os.path.isfile(os.path.join(path_one_sim,'ginette')):
    print ("ginette exist")
else:
    print ("ginette not exist")
    print("ginette will be compile in GINETTE_SENSI")
    print(" gfortran -o GINETTE_SENSI/ginette ../../src/ginette_V2.f")
    subprocess.call(["gfortran","-o","GINETTE_SENSI/ginette","../../src/ginette_V2.f"])  #creat

########### Setup
# /Users/mbp/Documents/my-project/python-snippets/notebook
f_param_bck=open("GINETTE_SENSI/E_parametre_backup.dat", "r")
setup_model=f_param_bck.read()
#time step in s
dt=900
#duration of the simulation in days
nb_day=2

#in meter
z_bottom=-0.40
z_top=0
dz=0.01



#Observation in meter
Obs1=-0.1
Obs2=-0.2
Obs3=-0.30
Obs4=-0.40

nb_zone=1
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



## write the parameters
nb_cell=-z_bottom/dz
cell1=-Obs1/dz
cell2=-Obs2/dz
cell3=-Obs3/dz
cell4=-Obs4/dz
state=0

setup_model=setup_model.replace('[dt]','%06.0fD+00' % dt)
setup_model=setup_model.replace('[nb_day]','%06.0f' % nb_day)
setup_model=setup_model.replace('[z_bottom]','%6.2e' % z_bottom)
setup_model=setup_model.replace('[az]','%7.3e' % -z_bottom)
setup_model=setup_model.replace('[dz]','%6.2e' % dz)
setup_model=setup_model.replace('[nb_cell]','%05.0f' % nb_cell)
setup_model=setup_model.replace('[cell1]','%05d' % cell1)
setup_model=setup_model.replace('[cell2]','%05d' % cell2)
setup_model=setup_model.replace('[cell3]','%05d' % cell3)
setup_model=setup_model.replace('[cell4]','%05d' % cell4)

setup_model=setup_model.replace('[state]','%1i' % state)

########### Zone
f_coor=open("GINETTE_SENSI/E_coordonnee.dat", "r")
f_zone=open("GINETTE_SENSI/E_zone.dat", 'w')
f_paramZ_bck=open("GINETTE_SENSI/E_zone_parameter_backup.dat", "r")
f_cdi_bck = open("GINETTE_SENSI/E_cdt_initiale_bck.dat", 'r')
f_cdl_bck = open("GINETTE_SENSI/E_cdt_aux_limites_bck.dat",'r')
f_paramZ_new = open("GINETTE_SENSI/E_zone_parameter.dat", 'w')
f_param_new = open("GINETTE_SENSI/E_parametre.dat", 'w')
f_cdi_new = open("GINETTE_SENSI/E_cdt_initiale.dat", 'w')
f_cdl_new = open("GINETTE_SENSI/E_cdt_aux_limites.dat",'w')
param_zone=f_paramZ_bck.read()
cdi=f_cdi_bck.read()
cdl=f_cdl_bck.read()
coord=pd.DataFrame()
coord = pd.read_csv(f_coor, names=["id", "x", "z"], header=None, delim_whitespace=True)



param_zone=param_zone.replace('[k1]','%8.2e' % val_k)
param_zone=param_zone.replace('[n1]','%6.2f' % val_n)
param_zone=param_zone.replace('[l1]','%6.2f' % val_l)
param_zone=param_zone.replace('[r1]','%6.2f' % val_r)



param_zone=param_zone.replace('[k2]','%8.2e' % val_k2)
param_zone=param_zone.replace('[r3]','%6.2f' % val_r3)






coord['zone'] =1

if nb_zone >= 2:
    coord['zone'] = np.where(coord['z'] <= thk2, 2,coord['zone'])
    coord['zone'] = np.where(coord['z'] <= thk3, 3,coord['zone'])

#display(coord)

f_param_new.write(setup_model)
f_paramZ_new.write(param_zone)
coord.zone.to_csv(f_zone, index = False, header=False)
f_zone.close()
f_paramZ_new.close()
f_paramZ_bck.close()

f_coor.close()
os.chdir(path_one_sim)
print("Ginette's steady")

# change parameters in input files for steady state
setup_model=setup_model.replace('[state]','%1i' % state)
ichi2=0
itempi=0
iclchgt=0
cdi=cdi.replace('[ichi2]','%1i' % ichi2)
cdi=cdi.replace('[itempi]','%1i' % itempi)
cdl=cdl.replace('[iclchgt]','%1i' % iclchgt)
f_cdl_new.write(cdl)
f_cdi_new.write(cdi)

f_param_bck.close()
f_param_new.close()
f_cdi_bck.close()
f_cdi_new.close()
f_cdl_bck.close()
f_cdl_new.close()

subprocess.call(["./ginette"])   


print("Ginette transient")
f_param_bck2=open("E_parametre_backup.dat", "r")
f_cdi_bck2 = open("E_cdt_initiale_bck.dat", 'r')
f_cdl_bck2 = open("E_cdt_aux_limites_bck.dat",'r')
f_param_new2 = open("E_parametre.dat", 'w')
f_cdi_new2 = open("E_cdt_initiale.dat", 'w')
f_cdl_new2 = open("E_cdt_aux_limites.dat",'w')
cdi=f_cdi_bck2.read()
cdl=f_cdl_bck2.read()
setup_model=f_param_bck2.read()
cdi=f_cdi_bck2.read()
cdl=f_cdl_bck2.read()
state=1
ichi2=1
itempi=1
iclchgt=1
cdi=cdi.replace('[ichi2]','%1i' % ichi2)
cdi=cdi.replace('[itempi]','%1i' % itempi)
cdl=cdl.replace('[iclchgt]','%1i' % iclchgt)
setup_model=setup_model.replace('[state]','%1i' % state)
f_cdl_new2.write(cdl)
f_cdi_new2.write(cdi)
f_param_new2.write(setup_model)

# paste pressure results from steady state
shutil.copyfile('S_pression_charge_temperature.dat','E_pression_initiale.dat')

f_param_bck2.close()
f_param_new2.close()
f_cdi_bck2.close()
f_cdi_new2.close()
f_cdl_bck2.close()
f_cdl_new2.close()



os.chdir(path_mini_lomos)# use copyfile()
