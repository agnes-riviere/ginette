#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 11:40:10 2024

@author: agnes_riviere antoine maison
"""
import os
import shutil
import subprocess
import numpy as np
import pandas as pd
import time
import threading
import matplotlib.pyplot as plt
import Functions_mesh_metis as fm


#%%
def convertMsh2Mail(meshName):    
    """
    @author agnes_riviere 

    Parameters
    ----------
    meshName : str
        Name of the mesh to be used by METIS.

    Returns
    -------
    None.

    """
    
    # libs_gfortran = ['gfortran']

    if os.path.isfile('mshTomail'):
        subprocess.call(["gfortran","-o","mshTomail","./mail_gmsh_2_met.f90"]) 
    else:
        print('compilation de mail_gmsh_2_met') 
        subprocess.call(["gfortran","-o","mshTomail","./mail_gmsh_2_met.f90"])     
        
    subprocess.call(["./mshTomail",meshName+".msh",meshName+".mail",meshName+".geom","2D/3D"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)    
    
    shutil.move(meshName+".mail","METIS/"+meshName+".mail")
    shutil.move("wl.n","METIS/wl.n")
    shutil.move("pm.n","METIS/pm.n")
    shutil.move("cable_nodes.dat","METIS/cable_nodes.dat")
    return
    

def editCOMM_METIS(meshName, df, n_iter):
    """
    

    Parameters
    ----------
    meshName : str
        Name of the mesh to be used by METIS.
    df : TYPE
        DESCRIPTION.
    n_iter : int
        Iteration number.

    Returns
    -------
    None.

    """
    with open("METIS/"+meshName+"_bck2.COMM","r") as f:
        param_model = f.read()
    
    param_model = param_model.replace('[Tskin]','%8.3f' % df.Tskin_CA[n_iter])
    
    with open("METIS/"+meshName+".COMM","w") as f:
        f.write(param_model)
    return


def preprocessMETIS(meshName, dicMETIS, caselabel="tempResults", visual=False):
    """
    @author agnes_riviere
    """
    respath = "Results/"+caselabel+"/"
    os.chdir("METIS/")
    
    # Data_file = False
    
    if meshName=='plat':
        fm.meshplat(meshName, L=dicMETIS["L"], xr1=dicMETIS["xr1"], xr2=dicMETIS["xr2"], z_bot=dicMETIS["z_bot"], z_sea=dicMETIS["z_sea"], z_sed=dicMETIS["z_sed"], 
                    cable_x=dicMETIS["cable_x"], cable_y=dicMETIS["Lb"], cable_r=dicMETIS["De"]/2, resolution_max=dicMETIS["resolution_max"], resolution_cable=dicMETIS["resolution_cable"], 
                    resolution_seabed=dicMETIS["resolution_seabed"], resolution_bound=dicMETIS["resolution_bound"], Dist_bound=dicMETIS["Dist_bound"],
                    Dist_seabed=dicMETIS["Dist_seabed"], Dist_cable=dicMETIS["Dist_cable"],visual=visual)
        fm.convertMsh2Mail(meshName)
        fm.extract_nodes_cable(meshName, "cable", dicMETIS["cable_x"], dicMETIS["Lb"], dicMETIS["De"]/2)
        node_in_line = fm.extract_nodes_lin_H_seabed(meshName, 'line', dicMETIS["z_sed"], dicMETIS["head_left"], dicMETIS["head_right"], dicMETIS["L"])
        node_in_line.head()
        # pl.plot_head(node_in_line)
        # v_left = extract_nodes_lin_H_V(meshName=meshName,bound_name="left_trans",x_bound=dicMETIS["xmin"], head_top=12.80001,head_bot=15,z_apply=2)
        # pl.plot_head(v_left)
        
        #fm.extract_nodes_lin_H_vert(meshName, "H_right", dicMETIS["L"], dicMETIS["head_right"], dicMETIS["head_bot"], dicMETIS["z_apply"]) 
        #fm.extract_nodes_lin_H_vert(meshName, "H_left", dicMETIS["L"], dicMETIS["head_top"], dicMETIS["head_left"], dicMETIS["z_apply"] )
        # Vérifiez que le fichier existe avant de l'ouvrir
        comm_file = meshName + "_bck.COMM"    
        if not os.path.isfile(comm_file):
            print(f"Erreur : Le fichier {comm_file} est introuvable.")
            sys.exit(1)  # Arrêtez le script avec un code d'erreur       
        fm.create_comm_metis_plat(meshName, dicMETIS)
        
        fm.extract_nodes_temp(meshName, dicMETIS["z_sed"], dicMETIS["T_sea"])
        
        os.chdir("../")
        return
    
    if meshName=='platsed':
        fm.meshsed(meshName, L=dicMETIS["L"], xr1=dicMETIS["xr1"], xr2=dicMETIS["xr2"], z_bot=dicMETIS["z_bot"], z_sed=dicMETIS["z_sed"],
                    cable_x=dicMETIS["cable_x"], cable_y=dicMETIS["Lb"], cable_r=dicMETIS["De"]/2, resolution_max=dicMETIS["resolution_max"], resolution_cable=dicMETIS["resolution_cable"], 
                   resolution_seabed=dicMETIS["resolution_seabed"], resolution_bound=dicMETIS["resolution_bound"], Dist_bound=dicMETIS["Dist_bound"],
                    Dist_seabed=dicMETIS["Dist_seabed"], Dist_cable=dicMETIS["Dist_cable"],visual=visual)
        fm.convertMsh2Mail(meshName)
        fm.extract_nodes_cable(meshName, "cable", dicMETIS["cable_x"], dicMETIS["Lb"], dicMETIS["De"]/2)
        fm.create_comm_metis_platsed(meshName, dicMETIS)
        
        os.chdir("../")
        
        return   
    
     
    if meshName=='platseddensitaire':
        fm.meshsed(meshName, L=dicMETIS["L"], xr1=dicMETIS["xr1"], xr2=dicMETIS["xr2"], z_bot=dicMETIS["z_bot"], z_sed=dicMETIS["z_sed"], 
                    cable_x=dicMETIS["cable_x"], cable_y=dicMETIS["Lb"], cable_r=dicMETIS["De"]/2, resolution_max=dicMETIS["resolution_max"], resolution_cable=dicMETIS["resolution_cable"], 
                   resolution_seabed=dicMETIS["resolution_seabed"], resolution_bound=dicMETIS["resolution_bound"], Dist_bound=dicMETIS["Dist_bound"],
                    Dist_seabed=dicMETIS["Dist_seabed"], Dist_cable=dicMETIS["Dist_cable"],visual=visual)
        fm.convertMsh2Mail(meshName)
        fm.extract_nodes_cable(meshName, "cable", dicMETIS["cable_x"], dicMETIS["Lb"], dicMETIS["De"]/2)
        fm.create_comm_metis_platseddens(meshName, dicMETIS)
        os.chdir("../")
        return     
    
    if meshName=='platdensitaire':
        fm.meshplat(meshName, L=dicMETIS["L"], xr1=dicMETIS["xr1"], xr2=dicMETIS["xr2"], z_bot=dicMETIS["z_bot"], z_sea=dicMETIS["z_sea"], z_sed=dicMETIS["z_sed"], 
                    cable_x=dicMETIS["cable_x"], cable_y=dicMETIS["Lb"], cable_r=dicMETIS["De"]/2, resolution_max=dicMETIS["resolution_max"], resolution_cable=dicMETIS["resolution_cable"], 
                    resolution_seabed=dicMETIS["resolution_seabed"], resolution_bound=dicMETIS["resolution_bound"], Dist_bound=dicMETIS["Dist_bound"],
                    Dist_seabed=dicMETIS["Dist_seabed"], Dist_cable=dicMETIS["Dist_cable"],visual=visual)
        fm.convertMsh2Mail(meshName)
        fm.extract_nodes_cable(meshName, "cable", dicMETIS["cable_x"], dicMETIS["Lb"], dicMETIS["De"]/2)
        node_in_line = fm.extract_nodes_lin_H_seabed(meshName, 'line', dicMETIS["z_sed"], dicMETIS["head_left"], dicMETIS["head_right"], dicMETIS["L"])
        node_in_line.head()
        # pl.plot_head(node_in_line)
        # v_left = extract_nodes_lin_H_V(meshName=meshName,bound_name="left_trans",x_bound=dicMETIS["xmin"], head_top=12.80001,head_bot=15,z_apply=2)
        # pl.plot_head(v_left)
        
        #fm.extract_nodes_lin_H_vert(meshName, "H_right", dicMETIS["L"], dicMETIS["head_right"], dicMETIS["head_bot"], dicMETIS["z_apply"]) 
        #fm.extract_nodes_lin_H_vert(meshName, "H_left", dicMETIS["L"], dicMETIS["head_top"], dicMETIS["head_left"], dicMETIS["z_apply"] )
        # Vérifiez que le fichier existe avant de l'ouvrir
        comm_file = meshName + "_bck.COMM"    
        if not os.path.isfile(comm_file):
            print(f"Erreur : Le fichier {comm_file} est introuvable.")
            sys.exit(1)  # Arrêtez le script avec un code d'erreur       
        fm.create_comm_metis_platdens(meshName, dicMETIS)
        
        fm.extract_nodes_temp(meshName, dicMETIS["z_sed"], dicMETIS["T_sea"])
        
        os.chdir("../")
        return
    
    
    print("Code runs through the second part of function preprocessMETIS")
    print("OK")
    print("#-----------------------#\n")
    return



def runAndReadMETIS(meshName, n_iter, dicMETIS):
    """
    Run METIS calculation and read the result files.

    Parameters
    ----------
    meshName : str
        Name of the mesh to be used by METIS.
    n_iter : int
        Iteration number.
    dicMETIS : dictionary
        Dictionary containing all the data of input file 'inputMETIS.dat', 
        describing the study case for METIS.

    Returns
    -------
    FLUX_METIS : float
        Heat flux through cable contour from METIS.
    TEMP_METIS : float
        Cable contour temperature from METIS.

    """
    print("METIS :")

    os.chdir("METIS/")



    
    os.system("./metis2020 < "+meshName+".COMM >/dev/null 2>&1")

    
    print("\nMETIS finished.")
    
    with open(meshName + ".tempprofc", 'r') as file:
        lines = file.readlines()

    # Trouver l'index de la dernière ligne commençant par #
    last_header_index = 0
    for i, line in enumerate(lines):
        if line.startswith('#'):
            last_header_index = i

    data = [
        line.strip().split()
        for line in lines[last_header_index + 1 :]
        if not line.startswith('#')
    ]

    # Créez un DataFrame à partir des données lues
    df_t = pd.DataFrame(data, columns=['co_x', 'co_y', 'temp'])


    # Assurez-vous que la colonne 'temp' est de type float
    df_t['temp'] = df_t['temp'].astype(float)

    # Calculer la moyenne de la colonne 'temp'
    TEMP_METIS = np.mean(df_t['temp'][:])
    
    print('TEMP_METIS:',TEMP_METIS)

    with open("flux_thermique.dat", 'r') as file:
        lines = file.readlines()

    # Trouver l'index de la dernière ligne commençant par #
    last_header_index = 0
    for i, line in enumerate(lines):
        if line.startswith('#'):
            last_header_index = i

    data = [
        line.strip().split()
        for line in lines[last_header_index + 1 :]
        if not line.startswith('#')
    ]
    # Créez un DataFrame à partir des données lues
    df_f = pd.DataFrame(data, columns=['co_x', 'co_y', 'flux'])

    # Assurez-vous que la colonne 'flux' est de type float
    df_f['flux'] = df_f['flux'].astype(float)

    # Calculer la somme de la colonne 'flux'
    FLUX_METIS = np.sum(df_f['flux'][:])
    print('FLUX_METIS:',FLUX_METIS)
    os.chdir("../")

    print("OK")
    print("#----------------------------#\n")
    return FLUX_METIS, TEMP_METIS

    
def moveMETISfiles(n_iter, meshName, caselabel='tempResults'):
    """
    Transport METIS results files from the calculation directory to the storage directory.

    Parameters
    ----------
    n_iter : int
        Iteration number.
    meshName : str
        Name of the mesh to be used by METIS.
    caselabel : str, optional
        Label of the current case, after which the storage folder is named. The default is 'tempResults'.

    Returns
    -------
    None.

    """
    METISdir = "Results/"+str(caselabel)+"/METIS/"
    
    shutil.move("METIS/"+meshName+".COMM", METISdir+meshName+"_"+str(n_iter)+".COMM")
    shutil.move("METIS/"+meshName+".htsoec.met", METISdir+meshName+".htsoec_"+str(n_iter)+".met")
    shutil.move("METIS/"+meshName+".tempprof", METISdir+meshName+"_"+str(n_iter)+".tempprof")
    shutil.move("METIS/"+meshName+".tempprofc", METISdir+meshName+"_"+str(n_iter)+".tempprofc")
    shutil.move("METIS/"+meshName+".vdaprof", METISdir+meshName+"_"+str(n_iter)+".vdaprof")
    # shutil.move("METIS/arete_flutxth.dat", "Results/METISResults/arete_flutxth_"+str(n_iter)+".dat")
    shutil.move("METIS/flux_thermique.dat", METISdir+"flux_thermique_"+str(n_iter)+".dat")
    shutil.move("METIS/"+meshName+".potprof", METISdir+meshName+"_"+str(n_iter)+".potprof")
    return
