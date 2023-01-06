#!/usr/bin/env python
import os
import sys
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import subprocess
libs_gfortran = ['gfortran']
# please compile ginette in the folder 1D_col
# path of the 1D_col directory
path_mini_lomos='/home/ariviere/Programmes/ginette/application/mini-LOMOS/'
os.chdir(path_mini_lomos)

path_one_sim='GINETTE_SENSI'

if os.path.isfile(os.path.join(path_one_sim,'ginette')):
    print ("ginette exist")
else:
    print ("ginette not exist")
    print("ginette will be compile in GINETTE_SENSI")
    print(" gfortran -o GINETTE_SENSI/ginette ../../src/ginette_V2.f")
    subprocess.call(["gfortran","-o","GINETTE_SENSI/ginette","../../src/ginette_V2.f"])  #creat

# Print the current working directory
print("Current working directory: {0}".format(os.getcwd()))

####################################################################### READ COMM FILE #######################################################################
File_com = "inversion.COMM"
Inversion_PT100 = "inversion_PT100.COMM"

if os.path.isfile(os.path.join(File_com)):
    print ("inversion.COMM exists")
else:
    print ("inversion.COMM not exist")
    print ("$1 : point name (no space) ex: Point1")
    print ("$2 : year (two last characters YY) ex for 2014: 14 ")
    print ("$3 : month (two last characters MM) ex for November: 11")
    print ("$4 : day (two last characters DD) ex: 07 ")
    print ("$5 : name temperature sensor (four first characters) ex: t502")
    print ("$6 : name pressure sensor (four first characters) ex: p502")
    print ("$7 : time step (s) ex: 900")
    print ("$8 : nb of observations ex: 3 if the 4 PT100 work, 2 if one is broken, etc. $8 = nb of working PT100-1; as the deepest working PT100 is used as a boundary condition. The upper boundary condition is the temperature sensor in the stream, in the pressure differential file.")
    print ("$9 : number of zones (clay, sand...) ex: 1 More than 2 would require more PT100 so it wouldn't make a lot of sense here.")
    print ("$10 : thickness of the upper zone (m) distance between the riverbed and the bottom of the upper zones. Type in 0 if there is just 1 zone.")
    print ("$11: maximum duration of the inversion (s) ex : 864000 (= 10 days)")
    
    
if os.path.isfile(os.path.join(Inversion_PT100)):
    print ("inversion_PT100.COMM exists")

else:
    print ("inversion_PT100.COMM not exist")
    print("each line each dz between 2 PT100")

with open(File_com, "r") as file_com:
    f_com = [i for line in file_com for i in line.split(' ') if i.strip()]
    COMM = f_com[0:12]
d = {ord(x):"" for x in "\n"}
COMM = [x.translate(d) for x in COMM]


print(COMM)


##################################### TREAT TME SERIES ##############################################################################################################################################
import subprocess
script_R=os.path.join(path_mini_lomos,'/formatToGinette.R')
form_file = subprocess.call("Rscript formatToGinette.R",shell=True)
                      
import numpy as np
import matplotlib.pyplot as plt

############################################################################ Generate parameters ########################################################################################################
inversion_parameter = "inversion_parameter.COMM"
if os.path.isfile(os.path.join(inversion_parameter)):
    print ("inversion_parameter.COMM exists")

else:
    print ("inversion_parameter.COMM not exist")
    print("Column 1 = parameters of interest: k = intrinsic permeability [m2], n = porosity , l = solid thermal conductivity [W m−1 K−1], r = solid density [kg m-3]")
    print("The bulk volumetric heat capacity of the porous medium is calculated  by the following equation : **c_mr_m = c_w r_w n + c_s r (1-n)**")
    print("c_w = specific heat capacity of water [334 000 J kg−1 K−1], r_w = water density [1 000 kg m-3], c_s = specific heat capacity of solid [J kg−1 K−1]")
    print("Column 2 = number of zones (e.g. clay, sand...)")
    print("Column 3 and 4 (respectively min and max) = range for the value test")
    print("Column 5 = number of tests within the previously defined range, Column 5 will be responsible for the number of simulations that will be launched.")
    print("-----------------------------------")
    print("example for 1 zone:")
    print("k 1 0001D-15 0001D-11 0003")
    print("n 1 0.40 0.40 0001")
    print("l 1 1300D-03 8400D-03 0002")
    print("r 1 2600D+00 2600D+00  0001")


import re

df_parameter = [] #as you want these as your headers 
with open(inversion_parameter) as parameter_com:
    for line in parameter_com:
        # remove whitespace at the start and the newline at the end
        line = line.strip()
        # split each column on whitespace
        col_param = re.split('\s+', line, maxsplit=5)
        df_parameter.append(col_param)

import pandas as pd
df = pd.DataFrame(df_parameter,columns=['parameter','zone','min','max','nb_values'])
df=df.set_index('parameter')

df = df.astype(float)

res_k=[]
res_n=[]
res_l=[]
res_r=[]
# Initializing variable
y = range(1,int(df.loc['k']['nb_values'])+1)
# Calculating result
if int(df.loc['k']['nb_values']) == 1:
    res_k.append(df.loc['k']['min'])
else : 
    res_k = np.logspace(np.log10(df.loc['k']['min']),np.log10(df.loc['k']['max']),int(df.loc['k']['nb_values']), endpoint = True)
if int(df.loc['n']['nb_values']) == 1:
    res_n.append(df.loc['n']['min'])
else :
    res_n = np.arange(df.loc['n']['min'],df.loc['n']['max'],int(df.loc['n']['nb_values']), endpoint = True)
if int(df.loc['l']['nb_values']) == 1:
    res_l.append(df.loc['l']['min'])
else :
    res_l = np.arange(df.loc['l']['min'],df.loc['l']['max'],int(df.loc['l']['nb_values']), endpoint = True)
if int(df.loc['r']['nb_values']) == 1:
    res_r.append(df.loc['r']['min'])
else :
    res_r = np.arange(df.loc['r']['min'],df.loc['r']['max'],int(df.loc['r']['nb_values']), endpoint = True)

# comment if you want run simulations
# Plotting the graph
""" fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(res_k, y, color = 'green')
ax.set_xscale("log")
plt.title('permeability')
plt.show() """
# porosity
""" y = range(1,int(df.loc['n']['nb_values'])+1)
fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(res_n, y, color = 'green')
plt.title('porosity')
plt.show() """
# thermal conductivity
""" y = range(1,int(df.loc['l']['nb_values'])+1)
fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(res_l, y, color = 'green')
plt.title('thermal conductivity')
plt.show() """
# heat capacity
""" y = range(1,int(df.loc['r']['nb_values'])+1)
fig, ax = plt.subplots(figsize = (9, 6))
ax.scatter(res_r, y, color = 'green')
plt.title('heat capacity')
plt.show() """

####### Parameter #####
print("Generating parameters' values")
lst_parameter=[]
from itertools import product
for x in product(res_k,res_n,res_l,res_r):
    lst_parameter.append(x)
np_parameter=np.array(lst_parameter)
print("your parameters are: k, n,l,r")
print(np_parameter)

ftestedvalues=os.path.join(path_one_sim,'tested_values')
with open(ftestedvalues, 'w') as f:
    f.writelines(' '.join(str(j) for j in i) + '\n' for i in lst_parameter)

####################################################################### position of PT100 #######################################################################
#Observation in meter
with open(Inversion_PT100, "r") as file1:
    f_list = [float(line) for line in file1 ]

Obs1=-f_list[0]/100
Obs2=Obs1-f_list[1]/100
Obs3=Obs2-f_list[2]/100
Obs4=Obs3-f_list[3]/100

#in meter
z_bottom=Obs4
z_top=0
dz=0.01

## write the parameters
nb_cell=-z_bottom/dz
cell1=-Obs1/dz
cell2=-Obs2/dz
cell3=-Obs3/dz
cell4=-Obs4/dz

#time step in s
dt=int(COMM[6])

#duration of the simulation in days
nb_day=int(COMM[10])/86400

nb_zone=int(COMM[8])

########### Setup
# /Users/mbp/Documents/my-project/python-snippets/notebook
f_param_bck=open("GINETTE_SENSI/E_parametre_backup.dat", "r")
setup_model=f_param_bck.read()

#####################################
## INITIALIZE PARAMETERS AND FILES ##
#####################################
# store name of file where parameters are stored
os.chdir(path_one_sim)
myModules = os.path.join(path_mini_lomos,path_one_sim)
if not myModules in sys.path:
    sys.path.append(myModules)
myModules
import  one_sim #
one_sim.HZ1D(dt,nb_day,z_bottom,dz,nb_cell,cell1,cell2,cell3,cell4,np_parameter)
os.chdir(path_mini_lomos)

script_R=os.path.join(path_mini_lomos,'Comparaison_mailles_sim-obs.R')
form_file = subprocess.call("Rscript Comparaison_mailles_sim-obs.R",shell=True)

