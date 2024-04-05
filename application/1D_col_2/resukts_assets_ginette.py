#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 10:44:01 2023

@author: th
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import datetime
import seaborn as sns
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

#%% Recharger data
os.chdir('/home/th/Documents/code/ginette/application/1D_col/sortie/data_inv_2/')
f1_c=np.loadtxt('c.txt')
f1_l=np.loadtxt('l.txt')
f1_a=np.loadtxt('a.txt')
f1_n=np.loadtxt('n.txt')
f1_s=np.loadtxt('s.txt')
f1_WT=np.loadtxt('WT.txt')
# f_RMSE=np.loadtxt('RMSE.txt')
f1_T5=np.loadtxt('T5.txt')
f1_T10=np.loadtxt('T10.txt')
f1_T30=np.loadtxt('T30.txt')
f1_T50=np.loadtxt('T50.txt')
f1_T70=np.loadtxt('T70.txt')
f1_T100=np.loadtxt('T100.txt')
f1_T150=np.loadtxt('T150.txt')
#%% Recharger data
os.chdir('/home/th/Documents/code/ginette/application/1D_col/sortie/data_inv_3/')
f2_c=np.loadtxt('c.txt')
f2_l=np.loadtxt('l.txt')
f2_a=np.loadtxt('a.txt')
f2_n=np.loadtxt('n.txt')
f2_s=np.loadtxt('s.txt')
f2_WT=np.loadtxt('WT.txt')
# f_RMSE=np.loadtxt('RMSE.txt')
f2_T5=np.loadtxt('T5.txt')
f2_T10=np.loadtxt('T10.txt')
f2_T30=np.loadtxt('T30.txt')
f2_T50=np.loadtxt('T50.txt')
f2_T70=np.loadtxt('T70.txt')
f2_T100=np.loadtxt('T100.txt')
f2_T150=np.loadtxt('T150.txt')
#%% Recharger data
os.chdir('/home/th/Documents/code/ginette/application/1D_col/sortie/data_inv_4_conductivite_high')
f3_c=np.loadtxt('c.txt')
f3_l=np.loadtxt('l.txt')
f3_a=np.loadtxt('a.txt')
f3_n=np.loadtxt('n.txt')
f3_s=np.loadtxt('s.txt')
f3_WT=np.loadtxt('WT.txt')
# f_RMSE=np.loadtxt('RMSE.txt')
f3_T5=np.loadtxt('T5.txt')
f3_T10=np.loadtxt('T10.txt')
f3_T30=np.loadtxt('T30.txt')
f3_T50=np.loadtxt('T50.txt')
f3_T70=np.loadtxt('T70.txt')
f3_T100=np.loadtxt('T100.txt')
f3_T150=np.loadtxt('T150.txt')

#%% Recharger data
os.chdir('/home/th/Documents/code/ginette/application/1D_col/sortie/data_inv_4_bis')
f4_c=np.loadtxt('c.txt')
f4_l=np.loadtxt('l.txt')
f4_a=np.loadtxt('a.txt')
f4_n=np.loadtxt('n.txt')
f4_s=np.loadtxt('s.txt')
f4_WT=np.loadtxt('WT.txt')
# f_RMSE=np.loadtxt('RMSE.txt')
f4_T5=np.loadtxt('T5.txt')
f4_T10=np.loadtxt('T10.txt')
f4_T30=np.loadtxt('T30.txt')
f4_T50=np.loadtxt('T50.txt')
f4_T70=np.loadtxt('T70.txt')
f4_T100=np.loadtxt('T100.txt')
f4_T150=np.loadtxt('T150.txt')
#%% Concatener les deux salves de modèles
f_c=np.concatenate((f3_c,f4_c))
f_l=np.concatenate((f3_l,f4_l))
f_a=np.concatenate((f3_a,f4_a))
f_n=np.concatenate((f3_n,f4_n))
f_s=np.concatenate((f3_s,f4_s))
f_WT=np.concatenate((f3_WT,f4_WT))
# f_RMSE=np.loadtxt('RMSE.txt')
f_T5=np.concatenate((f3_T5,f4_T5))
f_T10=np.concatenate((f3_T10,f4_T10))
f_T30=np.concatenate((f3_T30,f4_T30))
f_T50=np.concatenate((f3_T50,f4_T50))
f_T70=np.concatenate((f3_T70,f4_T70))
f_T100=np.concatenate((f3_T100,f4_T100))
f_T150=np.concatenate((f3_T150,f4_T150))
#%%
# f_c=f3_c
# f_l=f3_l
# f_a=f3_a
# f_n=f3_n
# f_s=f3_s
# f_WT=f3_WT
# f_T5=f3_T5
# f_T10=f3_T10
# f_T30=f3_T30
# f_T50=f3_T50
# f_T70=f3_T70
# f_T100=f3_T100
# f_T150=f3_T150

RMSE=np.ones(len(f_T5))
for i in range(0,len(f_T5)):
    RMSE5=np.sqrt(np.sum((f_T5[i,:]-data.loc[ind_j1:ind_j2,'5'].values)**2)/len(data.loc[ind_j1:ind_j2,'5']))
    RMSE10=np.sqrt(np.sum((f_T10[i,:]-data.loc[ind_j1:ind_j2,'10'].values)**2)/len(data.loc[ind_j1:ind_j2,'5']))
    RMSE30=np.sqrt(np.sum((f_T30[i,:]-data.loc[ind_j1:ind_j2,'30'].values)**2)/len(data.loc[ind_j1:ind_j2,'5']))
    RMSE50=np.sqrt(np.sum((f_T50[i,:]-data.loc[ind_j1:ind_j2,'50'].values)**2)/len(data.loc[ind_j1:ind_j2,'5']))
    RMSE70=np.sqrt(np.sum((f_T70[i,:]-data.loc[ind_j1:ind_j2,'70'].values)**2)/len(data.loc[ind_j1:ind_j2,'5']))
    RMSE100=np.sqrt(np.sum((f_T100[i,:]-data.loc[ind_j1:ind_j2,'100'].values)**2)/len(data.loc[ind_j1:ind_j2,'5']))
    RMSE150=np.sqrt(np.sum((f_T150[i,:]-data.loc[ind_j1:ind_j2,'150'].values)**2)/len(data.loc[ind_j1:ind_j2,'5']))

    RMSE[i]=np.sqrt(RMSE5**2+RMSE10**2+RMSE30**2+RMSE50**2+RMSE70**2+RMSE100**2+RMSE150**2)
    #RMSE[i]=np.exp(RMSE5)*np.exp(RMSE10)*np.exp(RMSE30)*np.exp(RMSE50)*np.exp(RMSE70)*np.exp(RMSE100)*np.exp(RMSE150)
    
#%% Calcul KGE 
# f_T5=f1_T5
mean_O5=np.mean(data.loc[ind_j1:ind_j2,'5'].values)
O5=data.loc[ind_j1:ind_j2,'5'].values
num5=np.ones(f_T5.shape[0])
denom5=np.ones(f_T5.shape[0])
r5=np.ones(f_T5.shape[0])
alpha5=np.ones(f_T5.shape[0])
beta5=np.ones(f_T5.shape[0])

mean_O10=np.mean(data.loc[ind_j1:ind_j2,'10'].values)
O10=data.loc[ind_j1:ind_j2,'10'].values
num10=np.ones(f_T10.shape[0])
denom10=np.ones(f_T10.shape[0])
r10=np.ones(f_T10.shape[0])
alpha10=np.ones(f_T10.shape[0])
beta10=np.ones(f_T10.shape[0])

mean_O30=np.mean(data.loc[ind_j1:ind_j2,'30'].values)
O30=data.loc[ind_j1:ind_j2,'30'].values
num30=np.ones(f_T30.shape[0])
denom30=np.ones(f_T30.shape[0])
r30=np.ones(f_T30.shape[0])
alpha30=np.ones(f_T30.shape[0])
beta30=np.ones(f_T30.shape[0])

mean_O50=np.mean(data.loc[ind_j1:ind_j2,'50'].values)
O50=data.loc[ind_j1:ind_j2,'50'].values
num50=np.ones(f_T50.shape[0])
denom50=np.ones(f_T50.shape[0])
r50=np.ones(f_T50.shape[0])
alpha50=np.ones(f_T50.shape[0])
beta50=np.ones(f_T50.shape[0])

mean_O70=np.mean(data.loc[ind_j1:ind_j2,'70'].values)
O70=data.loc[ind_j1:ind_j2,'70'].values
num70=np.ones(f_T70.shape[0])
denom70=np.ones(f_T70.shape[0])
r70=np.ones(f_T70.shape[0])
alpha70=np.ones(f_T70.shape[0])
beta70=np.ones(f_T70.shape[0])

mean_O100=np.mean(data.loc[ind_j1:ind_j2,'100'].values)
O100=data.loc[ind_j1:ind_j2,'100'].values
num100=np.ones(f_T100.shape[0])
denom100=np.ones(f_T100.shape[0])
r100=np.ones(f_T100.shape[0])
alpha100=np.ones(f_T100.shape[0])
beta100=np.ones(f_T100.shape[0])

mean_O150=np.mean(data.loc[ind_j1:ind_j2,'150'].values)
O150=data.loc[ind_j1:ind_j2,'150'].values
num150=np.ones(f_T150.shape[0])
denom150=np.ones(f_T150.shape[0])
r150=np.ones(f_T150.shape[0])
alpha150=np.ones(f_T150.shape[0])
beta150=np.ones(f_T150.shape[0])

mean_O=[mean_O5,mean_O10,mean_O30,mean_O50,mean_O70,mean_O100,mean_O150]
O=[O5,O10,O30,O50,O70,O100,O150]
num=[num5,num10,num30,num50,num70,num100,num150]
denom=[denom5,denom10,denom30,denom50,denom70,denom100,denom150]
r=[]

for i in range(0,f_T5.shape[0]):
    num5[i]=np.sum((f_T5[i,:]-np.mean(f_T5[i,:]))*(O5-mean_O5)) #Calcul numérateur coefficient de corrélation
    denom5[i]=np.sqrt(np.sum((f_T5[i,:]-np.mean(f_T5[i,:]))**2)*np.sum((O5-mean_O5)**2)) #Calcul dénominateur coefficient de corrélation
    alpha5[i]=np.std(f_T5[i,:])/np.std(O5)
    beta5[i]=np.mean(f_T5[i,:])/mean_O5
r5=num5/denom5   

KGE=1-np.sqrt((r5-1)**2+(alpha5-1)**2+(beta5-1)**2)

#%%
results=pd.DataFrame()
results['c']=f_c
results['lambda']=f_l
results['alpha']=f_a
results['n']=f_n
results['s']=f_s
results['WT']=f_WT
results['RMSE']=RMSE

sns.pairplot(results,hue='RMSE',palette='viridis')
#%% best results
plt.close('all')
indices=np.argsort(RMSE)
indices_inv=indices[::-1] #inverse l'indexation pour les plots

nb_best=10

best_results=pd.DataFrame()
best_results['c']=f_c[indices_inv[-nb_best:]]
best_results['lambda']=f_l[indices_inv[-nb_best:]]
best_results['alpha']=f_a[indices_inv[-nb_best:]]
best_results['n']=f_n[indices_inv[-nb_best:]]
best_results['s']=f_s[indices_inv[-nb_best:]]
best_results['WT']=f_WT[indices_inv[-nb_best:]]
best_results['RMSE']=RMSE[indices_inv[-nb_best:]]

sns.pairplot(best_results,hue='RMSE',palette='viridis')

#%% calcul teneur en eau avec les meilleurs paramètres de VG
# plt.close('all')

tr=best_results['s']
alpha=best_results['alpha']
n=best_results['n']
h=np.arange(0,100,0.01)
col=best_results['RMSE']
m=1-1/n

cmap=plt.cm.viridis

figi,ax=plt.subplots(1,1)
for ii in np.arange(0,nb_best,1):
    # if alpha[ii]<2 and n[ii]>1:
    theta=tr[ii]*(0.3-tr[ii])/(1+(alpha[ii]*h)**n[ii])**m[ii]
    ax.plot(theta,h,c=cmap(ii/len(n)))
# ax.set_xlim([0,0.5])
# ax.set_ylim([0,10])
plt.xlabel('Teneur en eau %')
plt.ylabel('Succion h en cm')
#%% Tri des paramètres de VG
for ii in range(len(best_results)):
    if best_results.loc[ii,'alpha']>2 and best_results.loc[ii,'n']<1:
        best_results.loc[ii,:]=np.nan

sns.pairplot(best_results,hue='RMSE',palette='viridis')

  
#%% comparer synth et data pour les X meileurs rmse
plt.close('all')

figi,ax=plt.subplots(2,2,sharex=True)
indices=np.argsort(RMSE)

for ii in np.arange(0,10,1):
    ax[0,0].plot(data.loc[ind_j1:ind_j2,'Timestamp'],f_T5[indices[ii],:])
    ax[0,1].plot(data.loc[ind_j1:ind_j2,'Timestamp'],f_T10[indices[ii],:])
    ax[1,0].plot(data.loc[ind_j1:ind_j2,'Timestamp'],f_T70[indices[ii],:])
    ax[1,1].plot(data.loc[ind_j1:ind_j2,'Timestamp'],f_T100[indices[ii],:])

ax[0,0].plot(data.loc[ind_j1:ind_j2,'Timestamp'],data.loc[ind_j1:ind_j2,'5'].values,'k')
ax[0,1].plot(data.loc[ind_j1:ind_j2,'Timestamp'],data.loc[ind_j1:ind_j2,'10'].values,'k')
ax[1,0].plot(data.loc[ind_j1:ind_j2,'Timestamp'],data.loc[ind_j1:ind_j2,'70'].values,'k')
ax[1,1].plot(data.loc[ind_j1:ind_j2,'Timestamp'],data.loc[ind_j1:ind_j2,'100'].values,'k')

ax[0,0].set_title('5 cm')
ax[0,0].set_ylabel('Température °C')
ax[0,1].set_title('10 cm')
ax[0,1].set_ylabel('Température °C')
ax[1,0].set_title('70 cm')
ax[1,0].set_ylabel('Température °C')
ax[1,1].set_title('100 cm')
ax[1,1].set_ylabel('Température °C')
ax[1,0].set_xlabel('Temps (MM-DD HH)')
ax[1,1].set_xlabel('Temps (MM-DD HH)')