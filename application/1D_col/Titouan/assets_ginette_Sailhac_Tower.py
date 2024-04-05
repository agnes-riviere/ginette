#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:27:47 2022

@author: th
"""

import pandas as pd
import numpy as np
import scipy.interpolate as interp
import matplotlib.pyplot as plt
import datetime
import matplotlib.dates as mdates
import matplotlib as mp
import os

from function_ginette_Sailhac_tower import preprocess,ginettout
#%%
data=pd.read_csv('/home/th/Documents/data/ASSETS/Temperatures.csv',sep=';',skiprows=2,header=None,names=['Time','Year','month','day','hour','minute','date','datH','Air','2.5','5','10','20','30','50','70','100','150','180'])
# data=pd.read_csv('/home/el/Codes/Titouan/data/Temperatures.csv',sep=';',skiprows=2,header=None,names=['Time','Year','month','day','hour','minute','date','datH','Air','2.5','5','10','20','30','50','70','100','150','180'])
data['Timestamp']=pd.to_datetime(data['Time'],dayfirst=True)
data=data.loc[(data['Timestamp']>='2018-09-05 00:00:00')&(data['Timestamp']<='2019-07-30 00:00:00'),:]
data=data.reset_index(drop=True)

#%% Paramètres du modèle
j1='2019-04-21 00:00:00'
j2='2019-04-24 23:45:00'
ind_j1=data.loc[data['Timestamp']==j1].index
ind_j2=data.loc[data['Timestamp']==j2].index
ind_j1=ind_j1[0]
ind_j2=ind_j2[0]
nb_j=round(pd.to_timedelta(data.loc[ind_j2,'Timestamp']-data.loc[ind_j1,'Timestamp'])/datetime.timedelta(days=1),1)


dt=900
nbcouche=1
nb_cell=preprocess(nb_j,dt,nbcouche)
#%% CONDITIONS INITIALES
#Température 
Prof=np.array([-0.025,-0.05,-0.1,-0.2,-0.3,-0.5,-0.7,-1,-1.5,-1.8])
f=interp.interp1d(Prof,np.array(data.iloc[ind_j1,9:19]), fill_value='extrapolate', kind='cubic')
Tini=pd.DataFrame(f(np.arange(0,-1.8,-0.01)))
Tini.to_csv('E_temperature_initiale.dat',header=None,index=None)
#Pression initiale
########### WT=-water_level.loc[ind_j1,1]/100 #m
WT=-1.5
rho_eau=1000 #kg.m-3
g = 9.81 #m.s-2
pression=np.round(np.array((WT-np.arange(0,-1.8,-0.01))*rho_eau*g),1)
p=pd.DataFrame(pression)
p.to_csv('E_pression_initiale.dat',index=None,header=None)

#%% CONDITIONS AUX LIMITES
#Température haute
T1=data.loc[ind_j1:ind_j2,'2.5'] # température premier capteur (2.5cm de profondeur)
T_air=data.loc[ind_j1:ind_j2,'Air'] # température de l'air
#Température basse
T9=data.loc[ind_j1:ind_j2,'180'] #température du dernier capteur (180 cm de profondeur)
#Création fichier conditions limites température
Texp=pd.DataFrame()
tt=np.arange(0,nb_j*24*3600,dt)
Texp['T']=T1
Texp['T2']=T9
Texp.to_csv('E_temp_t.dat',header=None,index=None,sep='\t')

#Flux/charge haut
flux_sup=6e-7
#Flux/charge bas
cdt_limite=pd.DataFrame(np.ones([len(data.loc[ind_j1:ind_j2]),2]),columns=['haut','bas'])
cdt_limite['haut']=flux_sup
################### cdt_limite['bas']=-water_level.loc[ind_j1:ind_j2,1].values/100
cdt_limite['bas']=6e-7
cdt_limite.to_csv('E_ec_bc_t.dat',index=None,header=None,sep='\t')

#%% GENERER GRAND NOMBRE DE MODELES
N=20000
res=[]
res_param=pd.DataFrame()
res_T=np.zeros([N,len(tt),7])
res_T5=pd.DataFrame()
res_T10=pd.DataFrame()
res_T30=pd.DataFrame()
res_T50=pd.DataFrame()
res_T70=pd.DataFrame()
res_T100=pd.DataFrame()
res_T150=pd.DataFrame()

dz=0.01
z_bottom=-1.8
x=np.arange(0,nb_j,900/86400)
y=np.arange(-dz/2,z_bottom,-dz)

k=0

for i in range(0,N):
    if i%200==0:
        print(i)
    capacity=np.array(np.random.uniform(500,6000,1))
    conductivity=np.array(np.random.uniform(0.01,11,1)) #Conductivité du solide jusqu'à 11 A.Rivière mail 23023/2023 10h01
    alpha=np.random.uniform(0.1,10,1)#inverse de l'entree dair [m]
    n=np.random.uniform(1.5,5,1)
    s=np.random.uniform(0.02,0.3,1)
    
    WT=float(np.random.uniform(0,-10,1))
    
    res_param.loc[i,'capacity']=float(capacity)
    res_param.loc[i,'conductivity']=float(conductivity)
    res_param.loc[i,'alpha']=alpha
    res_param.loc[i,'n']=n
    res_param.loc[i,'s']=s
    res_param.loc[i,'WT']=WT
    
    
    rho_eau=1000 #kg.m-3
    g = 9.81 #m.s-2
    pression=np.round(np.array((WT-np.arange(0,-1.8,-0.01))*rho_eau*g),1)
    p=pd.DataFrame(pression)
    p.to_csv('E_pression_initiale.dat',index=None,header=None)
    
    cdt_limite=pd.DataFrame(np.ones([len(data.loc[ind_j1:ind_j2]),2]),columns=['haut','bas'])
    cdt_limite['haut']=flux_sup
    ################### cdt_limite['bas']=-water_level.loc[ind_j1:ind_j2,1].values/100
    cdt_limite['bas']=WT
    cdt_limite.to_csv('E_ec_bc_t.dat',index=None,header=None,sep='\t')
    
    ginettout(capacity,conductivity,flux_sup,WT,nbcouche,alpha,n,s)
    
    temp=pd.read_table('Sim_temperature_profil_t.dat',header=None,delim_whitespace=True)
    flux=pd.read_table('Sim_heat_flux_profil_t.dat',header=None,delim_whitespace=True)
    if len(temp[2]) == len(x)*len(y):
        t=np.array(temp[2]).reshape(len(x),len(y))
        t=t.T
        T5=t[4,:]
        T10=t[9,:]
        T30=t[29,:]
        T50=t[49,:]
        T70=t[69,:]
        T100=t[99,:]
        T150=t[149,:]

    res_T[k,:,0]=T5
    res_T[k,:,1]=T10
    res_T[k,:,2]=T30
    res_T[k,:,3]=T50
    res_T[k,:,4]=T70
    res_T[k,:,5]=T100
    res_T[k,:,6]=T150
    
    # res_T5.loc[0]=T5
    # res_T10.loc[:,i]=T10
    # res_T30.loc[:,i]=T30
    # res_T50.loc[:,i]=T50
    # res_T70.loc[:,i]=T70
    # res_T100.loc[:,i]=T100
    # res_T150.loc[:,i]=T150



        # res.append([float(capacity),float(conductivity),alpha,n,s,WT,T5,T10,T30,T50,T70,T100,T150])
        
    if k==9:
        path_PS='/home/el/Codes/Titouan/ginette/application/1D_col/res/'
        path_TH='/home/th/Documents/code/ginette/application/1D_col/res/'
        
        f=open(path_TH+'c.txt','a')
        np.savetxt(f,res_param['capacity'].values)
        f.close()
        f=open(path_TH+'l.txt','a')
        np.savetxt(f,res_param['conductivity'].values)
        f.close()
        f=open(path_TH+'a.txt','a')
        np.savetxt(f,res_param['alpha'].values)
        f.close()
        f=open(path_TH+'n.txt','a')
        np.savetxt(f,res_param['n'].values)
        f.close()
        f=open(path_TH+'s.txt','a')
        np.savetxt(f,res_param['s'].values)
        f.close()
        f=open(path_TH+'WT.txt','a')
        np.savetxt(f,res_param['WT'].values)
        f.close()
        
        f=open(path_TH+'T5.txt','a')
        np.savetxt(f,res_T[:k+1,:,0])
        f.close()
        
        f=open(path_TH+'T10.txt','a')
        np.savetxt(f,res_T[:k+1,:,1])
        f.close()
        
        f=open(path_TH+'T30.txt','a')
        np.savetxt(f,res_T[:k+1,:,2])
        f.close()
        
        f=open(path_TH+'T50.txt','a')
        np.savetxt(f,res_T[:k+1,:,3])
        f.close()
        
        f=open(path_TH+'T70.txt','a')
        np.savetxt(f,res_T[:k+1,:,4])
        f.close()
        
        f=open(path_TH+'T100.txt','a')
        np.savetxt(f,res_T[:k+1,:,5])
        f.close()
        
        f=open(path_TH+'T150.txt','a')
        np.savetxt(f,res_T[:k+1,:,6])
        f.close()
        
        res_param=pd.DataFrame()
        res_T=np.zeros([N,len(tt),7])
        res_T5=pd.DataFrame()
        res_T10=pd.DataFrame()
        res_T30=pd.DataFrame()
        res_T50=pd.DataFrame()
        res_T70=pd.DataFrame()
        res_T100=pd.DataFrame()
        res_T150=pd.DataFrame()

        k=-1
    k=k+1



# #%% Recharger data
# f_c=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/c.txt')
# f_l=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/l.txt')
# f_a=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/a.txt')
# f_n=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/n.txt')
# f_s=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/s.txt')
# f_WT=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/WT.txt')
# f_T5=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/T5.txt')
# f_T10=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/T10.txt')
# f_T30=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/T30.txt')
# f_T50=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/T50.txt')
# f_T70=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/T70.txt')
# f_T100=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/T100.txt')
# f_T150=np.loadtxt('/home/th/Documents/code/ginette/application/1D_col/test/T150.txt')

