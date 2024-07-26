#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import gmsh
import numpy as np
import shutil  
import os    
import subprocess
import gmsh
import sys
import math
import matplotlib.pyplot as plt
from scipy.interpolate import griddata



def domain_state_time(z_top,z_bottom,dz,dt,state,nb_day,unsat):
    az=abs(z_top-z_bottom)
    # number of cell
    # ## write the parameters
    nb_cell=abs(z_top-z_bottom)/dz
#-----------------------------------------------------------------
    ## write the setup of the moddeled domain
    f_param_bck=open("E_parametre_backup.dat", "r")
    f_param_new = open("E_parametre.dat", 'w')
    setup_model=f_param_bck.read()
    setup_model=setup_model.replace('[dt]','%06.0fD+00' % dt)
    setup_model=setup_model.replace('[state]','%1i' % state)
    setup_model=setup_model.replace('[nb_day]','%06.0f' % nb_day)
    setup_model=setup_model.replace('[z_top]', '%7.3e' % z_top)
    setup_model=setup_model.replace('[dz]','%6.2e' % dz)
    setup_model=setup_model.replace('[nb_cell]','%05.0f' % nb_cell)
    setup_model=setup_model.replace('[unsat]','%1i' % unsat)
    setup_model=setup_model.replace('[dt]','%06.0fD+00' % dt)
    setup_model=setup_model.replace('[z_bottom]', '%7.3e' % z_bottom)
    setup_model=setup_model.replace('[az]','%7.3e' % az)



    f_param_new.write(setup_model)
    f_param_bck.close()
    f_param_new.close()


    return nb_cell



def geometry(z_top,z_bottom,dz,nb_zone):    
        ########### Zone of parameters
    f_coor=open("E_coordonnee.dat", "w")
    f_zone=open("E_zone.dat", 'w')
    coord=pd.DataFrame()
    # calculate the coordinates to know the number of cell in the domain    
    # Coordinate  
    zvalues =  np.sort(np.arange(z_bottom+dz/2,z_top,dz ))[::-1]

    xvalues = np.array([0.5]) # 1D column
    zz, xx = np.meshgrid(zvalues, xvalues)
    NT = np.prod(zz.shape)
    data = {
        "x": np.reshape(xx,NT),
        "z": np.reshape(zz,NT)}
    coord = pd.DataFrame(data=data)
    coord['id']=coord.index.values.astype(int)
    coord['id']=coord['id']+1
    cols = coord.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    coord = coord[cols] 
    coord.to_csv(f_coor, index = False, sep=' ', header=False)
    #zone parameter by cell ((homogenous domain = 1 zone))
    coord['zone'] =1
    #Pour plusieurs zones modification TH
    # if nb_zone >= 2:
    #     for i in range(2,int(nb_zone)+1):
    #         coord['zone'] = np.where(coord['z'] <= coord.loc[((i-1)*z_top*100/nb_zone),'z'], i,coord['zone'])
            
    
    #coord['zone'] = np.where(coord['z'] <= coord.loc[i-1,'z'], i,coord['zone'])
            #coord['zone'] = np.where(coord['z'] <= thk2, 2,coord['zone'])



    coord.zone.to_csv(f_zone, index=False, header=False)
    
    # close files    
    f_zone.close()
    f_coor.close()
    return coord


def initial_boundary_condition(val_WT,itlecture):
    # Initial conditions
    f_IC_bck=open("E_cdt_initiale_backup.dat","r")
    IC_model=f_IC_bck.read()
    IC_model=IC_model.replace('[head_ini]', '%05.2fD+00' % val_WT)

    ## write the boundary conditions
    f_bc_bck=open("E_cdt_aux_limites_backup.dat", "r")
    bc_model=f_bc_bck.read()
    bc_model=bc_model.replace('[top]', '%08.2fD+00' % val_WT)
    bc_model=bc_model.replace('[itlecture]', '%8i' % itlecture)   
    f_bc_new = open("E_cdt_aux_limites.dat", 'w')
    f_IC_new=open("E_cdt_initiale.dat","w")
    #Write new ginette files
    f_IC_new.write(IC_model)
    f_bc_new.write(bc_model)
    f_IC_new.close()
    f_IC_bck.close()
    f_bc_bck.close()
    f_bc_new.close()
    
def parameter_zone_geometrie(nb_zone,def_zone,dz,z_top,z_bottom):
    f_zone=open("E_zone.dat", 'w')
            ########### Zone of parameters
    f_coor=open("E_coordonnee.dat", "w")
    coord=pd.DataFrame()
    # calculate the coordinates to know the number of cell in the domain   
    if z_top>z_bottom:
        zvalues = np.arange(z_top-dz/2, z_bottom,-dz );
    else:
        zvalues = np.arange(z_bottom-dz/2, z_top,-dz );    
        
    xvalues = np.array([0.5]);
    zz, xx = np.meshgrid(zvalues, xvalues)
    NT = np.prod(zz.shape)
    data = {
        "x": np.reshape(xx,NT),
        "z": np.reshape(zz,NT)}
    coord = pd.DataFrame(data=data)
    coord['id']=coord.index.values.astype(int)
    coord['id']=coord['id']+1
    cols = coord.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    coord = coord[cols] 
    coord.to_csv(f_coor, index = False, sep=' ', header=False)
    #zone parameter by cell ((homogenous domain = 1 zone))
    coord['zone'] =1
    #Pour plusieurs zones modification 
    # def_zone bottom of each zone if z below def_zone[i] then zone = i+1
    if (nb_zone >= 2):
        for i in range(2,int(nb_zone)+1):
            coord['zone'] = np.where(coord['z'] <= def_zone[i-1], i,coord['zone'])

    #Write new ginette files
    coord.zone.to_csv(f_zone, index = False, header=False)
    
    # close files    
    f_zone.close() 
   
    
def parameter_zone_value(zone, val_k, val_n, val_a, val_nVG, val_swres):
    # Chemin du fichier où écrire les paramètres
    fichier_parametres = "E_zone_parameter.dat"
    # Construction de la ligne de paramètres pour la zone spécifiée
    ligne_parametres = f'{zone} {val_k:8.2e} {val_n:6.2f} {val_a:6.2f} {val_nVG:6.2f} {val_swres:6.2f}\n'
    
    # Écriture de la ligne de paramètres à la fin du fichier
    with open(fichier_parametres, 'a') as f:
        f.write(ligne_parametres)