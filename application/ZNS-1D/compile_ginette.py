#!/usr/bin/env python

import os
import numpy as np
from pathlib import Path
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib as mpl
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
dz=0.0035

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
val_a=2.70000 #m-1 alpha_vg
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
# Initial conditions

f_IC_bck=open("E_cdt_initiale_backup.dat","r")

IC_model=f_IC_bck.read()
IC_model=IC_model.replace('[head_ini]', '%05.0fD+00' % top)
## write the boundary conditions
f_bc_bck=open("E_cdt_aux_limites_backup.dat", "r")
bc_model=f_bc_bck.read()

bc_model=bc_model.replace('[top]', '%08.0fD+00' % top)
bc_model=bc_model.replace('[bot]','%08.0fD+00' % bot)

########### Zone
f_coor=open("E_coordonnee.dat", "w")
f_zone=open("E_zone.dat", 'w')
f_paramZ_bck=open("E_zone_parameter_backup.dat", "r")
f_paramZ_new = open("E_zone_parameter.dat", 'w')
f_param_new = open("E_parametre.dat", 'w')
f_bc_new = open("E_cdt_aux_limites.dat", 'w')
f_IC_new=open("E_cdt_initiale.dat","w")
param_zone=f_paramZ_bck.read()
coord=pd.DataFrame()    
#coord = pd.read_csv(f_coor, names=["id", "x", "z"], header=None, delim_whitespace=True)



param_zone=param_zone.replace('[k1]','%8.2e' % val_k)
param_zone=param_zone.replace('[n1]','%6.2f' % val_n)
param_zone=param_zone.replace('[r1]','%6.2f' % val_r)
param_zone=param_zone.replace('[a1]','%8.2e' % val_a)
param_zone=param_zone.replace('[nVG1]','%6.2f' % val_nVG)
param_zone=param_zone.replace('[swres1]','%6.2f' % val_swres)



zvalues = np.arange(dz/2, z_top,dz );
xvalues = np.array([0.5]);
zz, xx = np.meshgrid(zvalues, xvalues)
NT = np.product(zz.shape)

data = {
    "x": np.reshape(xx,NT),
    "z": np.reshape(zz,NT)
}
coord = pd.DataFrame(data=data)
coord['id']=coord.index.values.astype(int)
coord['id']=coord['id']+1
cols = coord.columns.tolist()
cols = cols[-1:] + cols[:-1]
coord = coord[cols] 
coord.to_csv(f_coor, index = False, sep=' ', header=False)

#zone parameter by cell
coord['zone'] =1


#display(coord)
f_IC_new.write(IC_model)
f_param_new.write(setup_model)
f_paramZ_new.write(param_zone)
f_bc_new.write(bc_model)
coord.zone.to_csv(f_zone, index = False, header=False)
f_zone.close()
f_coor.close()
f_paramZ_new.close()
f_paramZ_bck.close()
f_param_bck.close()
f_param_new.close()
f_IC_new.close()
f_bc_bck.close()
f_bc_new.close()
f_coor.close()

subprocess.call(["./ginette"])   

saturation_profile = pd.read_table('S_saturation_profil_t.dat',delim_whitespace=True,header=None)
saturation_profile.columns=[ "time",  "z","sat"]
#print(saturation_profile.head())
saturation_profile_s = pd.read_table('/home/ariviere/Documents/Encadrements/2021_These_Ramon/data/data_fig_4_5_6/sandyclay_WT25_prop_model2.txt',delim_whitespace=True,header=None,skiprows=[0,1,2])
saturation_profile_s.columns=[ "z",  "sat","rho","vp","vs"]
saturation_profile_s.z=saturation_profile_s.z+40
plt.figure()
plt.style.use('seaborn')

plt.scatter(saturation_profile.sat, saturation_profile.z, s=10, alpha=1, color='mediumblue',marker='.')
plt.scatter(saturation_profile_s.sat,saturation_profile_s.z, s=5, c='r', marker=",")
plt.xlabel('Saturation')
plt.ylabel('z (m)')
plt.show()




pressure_profile = pd.read_table('S_pressure_profil_t.dat',delim_whitespace=True,header=None)
pressure_profile.columns=[ "time",  "z","p","h"]
print(saturation_profile_s)

#plt.figure()
#plt.style.use('seaborn')

#plt.scatter(pressure_profile.p,pressure_profile.z, s=10, alpha=1, color='mediumblue',marker='.')
#plt.xlabel('Pressure (Pa)')
#plt.ylabel('z (m)')
#plt.show()