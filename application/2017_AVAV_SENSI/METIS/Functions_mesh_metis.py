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
import matplotlib.colors as colors
from scipy.interpolate import griddata

#%%
# Pour chaque élément, tracez une ligne entre chaque paire de nœuds
def get_coord(noeud,n):
    series = noeud[noeud['identifiant'] == n]
    if not series.empty:
        return series['x'].values[0], series['y'].values[0]
    else:
        return None, None


def createInputDictionary(inputFile):
    """
    Parameters
    ----------
    inputFile : ASCII text file
        ASCII text file containing METIS input data

    Returns
    -------
    dictionary : 
        Creation of a dictionary based on input text file data
    """
         
    df  =  pd.read_csv(inputFile, delimiter = ' = ',header = 0, comment = '#',names = ['data', 'value'],engine='python')
    
    dictionary  =  df.set_index('data')['value'].to_dict()
        
    return dictionary

def find_file(filename, search_path):
    for root, dir, files in os.walk(search_path):
        if filename in files:
            return os.path.join(root, filename)
    raise FileNotFoundError(f"Fichier {filename} non trouvé dans {search_path}")

def run_metis(meshName):    
    """

    Parameters
    ----------
    meshName : mesh and job name for the METIS computation

    Returns
    -------
    The "meshName".msh file 
    
    Conversion of gmsh mesh and METIS mesh

    """
    # ver
    # libs_gfortran = ['gfortran']

    if os.path.isfile('mshTomail'):
        # print("convertisseur msh to mail")
        subprocess.call(["gfortran","-o","mshTomail","./mail_gmsh_2_met.f90"]) 
        
    else:
        # print('compilation de mail_gmsh_2_met') 
        subprocess.call(["gfortran","-o","mshTomail","./mail_gmsh_2_met.f90"])     
        
    subprocess.call(["./mshTomail",meshName+".msh",meshName+".mail",meshName+".geom","2D/3D"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)    
    return


def create_comm_metis_platsed(meshName,dicMETIS):
    job_name = meshName
	

    z_sed = dicMETIS["z_sed"]
    ywater = (dicMETIS["z_sea"]-dicMETIS["z_sed"])*0.5
    z_botmean = 0.5*(dicMETIS["z_bot"]+dicMETIS["z_sed"])
    xmean = 0.5*(dicMETIS["xmin"]+dicMETIS["L"])	
    #### ouverture du fichier comm backup
    # f_param_bck = open("plat_bck.COMM","r")
    
    f_param_bck = open(meshName+"_bck.COMM","r")
    
    #### lecture du fichier comm backup
    param_model = f_param_bck.read()
    #### modify commande file METIS
    param_model = param_model.replace('[job_name]', job_name)
    param_model = param_model.replace('[Tsea]', '%8.2f' % dicMETIS["T_sea"])
    param_model = param_model.replace('[Tbot]', '%8.2f' % dicMETIS["T_bot"])
    param_model = param_model.replace('[akh1]','%8.2e' % dicMETIS["val_Kh1"])
    param_model = param_model.replace('[akv1]','%8.2e' % dicMETIS["val_Kv1"])
    param_model = param_model.replace('[l1]','%8.2e' % dicMETIS["lambda_pm"])
    param_model = param_model.replace('[xmin_geom]','%8.5f' % dicMETIS["xmin_geom"])
    param_model = param_model.replace('[xmax_geom]','%8.5f' % (dicMETIS["L"]+1.))
    param_model = param_model.replace('[ymin_geom]','%8.5f' % dicMETIS["ymin_geom"])
    param_model = param_model.replace('[ymax_geom]','%8.5f' % dicMETIS["ymax_geom"])
    param_model = param_model.replace('[z_sed2]','%8.2f' % dicMETIS['z_sed'])
    param_model = param_model.replace('[xmax]','%8.2f' % dicMETIS["L"])
    param_model = param_model.replace('[xmin]','%8.2f' % dicMETIS["xmin"])
    param_model = param_model.replace('[z_sed1]','%8.2f' % z_sed)
    param_model = param_model.replace('[z_bot]','%8.2f' % dicMETIS["z_bot"])
    param_model = param_model.replace('[xmean]','%8.2f' % xmean)
    param_model = param_model.replace('[z_botmean]','%8.2f' % z_botmean)
    param_model = param_model.replace('[z_sed]','%8.2f' % dicMETIS["z_sed"])

    f_param_bck.close()

      # Imprimer le répertoire courant
    print("Current directory:", os.getcwd())
    print("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
    with open(job_name+"_bck2.COMM","w") as f:
        f.write(param_model)


    return



def create_comm_metis_platseddens(meshName,dicMETIS):
    job_name = meshName
    z_sed = dicMETIS["z_sed"]
    ywater = (dicMETIS["z_sea"]-dicMETIS["z_sed"])*0.5
    z_botmean = 0.5*(dicMETIS["z_bot"]+dicMETIS["z_sed"])
    xmean = 0.5*(dicMETIS["xmin"]+dicMETIS["L"])	
    #### ouverture du fichier comm backup
    # f_param_bck = open("plat_bck.COMM","r")
    
    f_param_bck = open(meshName+"_bck.COMM","r")
    
    #### lecture du fichier comm backup
    param_model = f_param_bck.read()
    #### modify commande file METIS
    param_model = param_model.replace('[job_name]', job_name)
    param_model = param_model.replace('[hsea]', '%8.2f' % dicMETIS["head_top"])
    param_model = param_model.replace('[Tsea]', '%8.2f' % dicMETIS["T_sea"])
    param_model = param_model.replace('[Tbot]', '%8.2f' % dicMETIS["T_bot"])
    param_model = param_model.replace('[akh1]','%8.2e' % dicMETIS["val_Kh1"])
    param_model = param_model.replace('[akv1]','%8.2e' % dicMETIS["val_Kv1"])
    param_model = param_model.replace('[akh2]','%8.2e' % dicMETIS["val_Kh2"])
    param_model = param_model.replace('[akv2]','%8.2e' % dicMETIS["val_Kv2"])
    param_model = param_model.replace('[l1]','%8.2e' % dicMETIS["lambda_pm"])
    param_model = param_model.replace('[xmin_geom]','%8.5f' % dicMETIS["xmin_geom"])
    param_model = param_model.replace('[xmax_geom]','%8.5f' % (dicMETIS["L"]+1.))
    param_model = param_model.replace('[ymin_geom]','%8.5f' % dicMETIS["ymin_geom"])
    param_model = param_model.replace('[ymax_geom]','%8.5f' % dicMETIS["ymax_geom"])
    param_model = param_model.replace('[z_sed]','%8.2f' % dicMETIS['z_sed'])
    param_model = param_model.replace('[xmax]','%8.2f' % dicMETIS["L"])
    param_model = param_model.replace('[xmin]','%8.2f' % dicMETIS["xmin"])
    param_model = param_model.replace('[z_sed1]','%8.2f' % z_sed)
    param_model = param_model.replace('[z_bot]','%8.2f' % dicMETIS["z_bot"])
    param_model = param_model.replace('[xmean]','%8.2f' % xmean)
    param_model = param_model.replace('[z_botmean]','%8.2f' % z_botmean)
    param_model = param_model.replace('[z_sed]','%8.2f' % dicMETIS["z_sed"])

    f_param_bck.close()

    with open(job_name+"_bck2.COMM","w") as f:
        f.write(param_model)

    return


def create_comm_metis_plat(meshName,dicMETIS):
    job_name = meshName
    yseabed1 = dicMETIS["z_sed"]
    yseabed2 = dicMETIS["z_sed"]
    ywater = (dicMETIS["z_sea"]-dicMETIS["z_sed"])*0.5
    z_botmean = 0.5*(dicMETIS["z_bot"]+dicMETIS["z_sed"])
    xmean = 0.5*(dicMETIS["xmin"]+dicMETIS["L"])
    # Vérifiez que le fichier existe avant de l'ouvrir
    comm_file = job_name + "_bck.COMM"    
    if not os.path.isfile(comm_file):
        print(f"Erreur : Le fichier {comm_file} est introuvable.")
        sys.exit(1)  # Arrêtez le script avec un code d'erreur
    
    with open(job_name + "_bck.COMM", "r") as f_param_bck:
        #test find file comm and if not error message and stop
        print("file comm found")

        param_model = f_param_bck.read()
        #### modify com file METIS
        param_model = param_model.replace('[job_name]', job_name)
        param_model = param_model.replace('[Tsea]', '%8.2f' % dicMETIS["T_sea"])
        param_model = param_model.replace('[Tbot]', '%8.2f' % dicMETIS["T_bot"])
        param_model = param_model.replace('[akh1]','%8.2e' % dicMETIS["val_Kh1"])
        param_model = param_model.replace('[akv1]','%8.2e' % dicMETIS["val_Kv1"])
        param_model = param_model.replace('[akh2]','%8.2e' % dicMETIS["val_Kh2"])
        param_model = param_model.replace('[akv2]','%8.2e' % dicMETIS["val_Kv2"])
        param_model = param_model.replace('[l1]','%8.2e' % dicMETIS["lambda_pm"])
        param_model = param_model.replace('[xmin_geom]','%8.5f' % dicMETIS["xmin_geom"])
        param_model = param_model.replace('[xmax_geom]','%8.5f' % (dicMETIS["L"]+1.))
        param_model = param_model.replace('[ymin_geom]','%8.5f' % dicMETIS["ymin_geom"])
        param_model = param_model.replace('[ymax_geom]','%8.5f' % dicMETIS["ymax_geom"])
        param_model = param_model.replace('[z_sea]','%8.2f' % dicMETIS["z_sea"])
        param_model = param_model.replace('[yseabed2]','%8.2f' % yseabed2)
        param_model = param_model.replace('[xmax]','%8.2f' % dicMETIS["L"])
        param_model = param_model.replace('[ywater]','%8.2f' % ywater)
        param_model = param_model.replace('[xmin]','%8.2f' % dicMETIS["xmin"])
        param_model = param_model.replace('[yseabed1]','%8.2f' % yseabed1)
        param_model = param_model.replace('[z_bot]','%8.2f' % dicMETIS["z_bot"])
        param_model = param_model.replace('[xmean]','%8.2f' % xmean)
        param_model = param_model.replace('[Qr]','%8.2e' % dicMETIS["Qr"])
        param_model = param_model.replace('[Ql]','%8.2e' % dicMETIS["Ql"])
        param_model = param_model.replace('[z_botmean]','%8.2f' % z_botmean)
        param_model = param_model.replace('[z_sed]','%8.2f' % dicMETIS["z_sed"])
        param_model = param_model.replace('[l2]','%8.2e' % dicMETIS["lambda_wl"])

    with open(job_name+"_bck2.COMM","w") as f:
        f.write(param_model)

    return



def create_comm_metis_platdens(meshName,dicMETIS):
    job_name = meshName
    print("job_name",job_name)
    yseabed1 = dicMETIS["z_sed"]
    yseabed2 = dicMETIS["z_sed"]
    ywater = (dicMETIS["z_sea"]-dicMETIS["z_sed"])*0.5
    z_botmean = 0.5*(dicMETIS["z_bot"]+dicMETIS["z_sed"])
    xmean = 0.5*(dicMETIS["xmin"]+dicMETIS["L"])
    # Vérifiez que le fichier existe avant de l'ouvrir
    comm_file = job_name + "_bck.COMM"    
    if not os.path.isfile(comm_file):
        print(f"Erreur : Le fichier {comm_file} est introuvable.")
        sys.exit(1)  # Arrêtez le script avec un code d'erreur
    
    with open(job_name + "_bck.COMM", "r") as f_param_bck:
        #test find file comm and if not error message and stop
        print("file comm found")

        param_model = f_param_bck.read()
        #### modify com file METIS
        param_model = param_model.replace('[job_name]', job_name)
        param_model = param_model.replace('[Tsea]', '%8.2f' % dicMETIS["T_sea"])
        param_model = param_model.replace('[Tbot]', '%8.2f' % dicMETIS["T_bot"])
        param_model = param_model.replace('[akh1]','%8.2e' % dicMETIS["val_Kh1"])
        param_model = param_model.replace('[akv1]','%8.2e' % dicMETIS["val_Kv1"])
        param_model = param_model.replace('[akh2]','%8.2e' % dicMETIS["val_Kh2"])
        param_model = param_model.replace('[akv2]','%8.2e' % dicMETIS["val_Kv2"])
        param_model = param_model.replace('[l1]','%8.2e' % dicMETIS["lambda_pm"])
        param_model = param_model.replace('[xmin_geom]','%8.5f' % dicMETIS["xmin_geom"])
        param_model = param_model.replace('[xmax_geom]','%8.5f' % (dicMETIS["L"]+1.))
        param_model = param_model.replace('[ymin_geom]','%8.5f' % dicMETIS["ymin_geom"])
        param_model = param_model.replace('[ymax_geom]','%8.5f' % dicMETIS["ymax_geom"])
        param_model = param_model.replace('[z_sea]','%8.2f' % dicMETIS["z_sea"])
        param_model = param_model.replace('[yseabed2]','%8.2f' % yseabed2)
        param_model = param_model.replace('[xmax]','%8.2f' % dicMETIS["L"])
        param_model = param_model.replace('[ywater]','%8.2f' % ywater)
        param_model = param_model.replace('[xmin]','%8.2f' % dicMETIS["xmin"])
        param_model = param_model.replace('[yseabed1]','%8.2f' % yseabed1)
        param_model = param_model.replace('[z_bot]','%8.2f' % dicMETIS["z_bot"])
        param_model = param_model.replace('[xmean]','%8.2f' % xmean)
        param_model = param_model.replace('[Qr]','%8.2e' % dicMETIS["Qr"])
        param_model = param_model.replace('[Ql]','%8.2e' % dicMETIS["Ql"])
        param_model = param_model.replace('[z_botmean]','%8.2f' % z_botmean)
        param_model = param_model.replace('[z_sed]','%8.2f' % dicMETIS["z_sed"])
        param_model = param_model.replace('[l2]','%8.2e' % dicMETIS["lambda_wl"])

    with open(job_name+"_bck2.COMM","w") as f:
        f.write(param_model)

    return




def extract_nodes_lin_H_hor(meshName, bound_name, visual,datFile, HdatFile):
    print(os.getcwd())
    filename = meshName + ".mail"
    try:
        file_path = find_file(filename, os.getcwd())
        print(f"Fichier {filename} trouvé dans {file_path}")
    except FileNotFoundError as e:
        print(e)
        return

    # Lire le fichier '.mail' et extraire les coordonnées
    with open(file_path, 'r') as file:
        content = file.read()
        lines = content.split('\n')
        set_up = lines[0]
        nb_node_elt = set_up.split() 

    nb_node = int(nb_node_elt[0])
    nb_elt = int(nb_node_elt[1])


    # Vérifier l'existence des fichiers
    if not os.path.exists(datFile):
        print(f"Le fichier {datFile} n'existe pas.")
        return
    if not os.path.exists(HdatFile):
        print(f"Le fichier {HdatFile} n'existe pas.")
        return
    # Lire datFile
    df_seabed = pd.read_csv(datFile, sep=r'\s+', names=['x', 'y'], engine='python')
    #order df_seabed by x   
    df_seabed.sort_values(by=['x'], inplace=True)

    xmin, xmax = df_seabed['x'].min(), df_seabed['x'].max()
    ymin, ymax = df_seabed['y'].min(), df_seabed['y'].max()    

    # Lire HdatFile
    df_freesurf = pd.read_csv(HdatFile, sep=r'\s+', names=['x', 'freesurf'],engine='python')       
    df_freesurf.sort_values(by=['x'], inplace=True)
    df_freesurf['H'] = df_freesurf['freesurf'] - df_seabed['y']    

    # Lire le fichier 'job1.mail' et extraire les coordonnées
    with open(filename, 'r') as file:
        content = file.read()
        lines = content.split('\n')
        node_lines = lines[1:nb_node+1]
        coordinates = [line.split() for line in node_lines]
        x = [float(coord[0]) for coord in coordinates]
        y = [float(coord[1]) for coord in coordinates]
        noeud = pd.DataFrame({'x': x, 'y': y})

    # Réinitialiser l'index et renommer la colonne
    noeud.reset_index(inplace=True)
    noeud.rename(columns={'index': 'identifiant'}, inplace=True)

    # Mettre à jour le dataframe 'noeud'
    noeud['identifiant'] += 1
    model_xmin, model_xmax = noeud['x'].min(), noeud['x'].max()
    print('Model=   xmin:', model_xmin, 'xmax:', model_xmax)

    # add the point x=0 and y=df_seabed['y'][0] and x=xmax and y=df_seabed['y'][-1]
    # create dataframe with two rows
    points_to_add = pd.DataFrame({'x': [model_xmin, model_xmax,xmax-0.01], 'y': [df_seabed['y'].iloc[0], df_seabed['y'].iloc[-1], df_seabed['y'].iloc[-1]]})
    # add the two rows to df_seabed
    df_seabed = pd.concat([df_seabed, points_to_add], ignore_index=True)
#   create a x betwen xmin and xmax each 0.01
    x = np.arange(model_xmin, model_xmax, 0.01)
     # interpolation spline   
    spline_seabed = pd.DataFrame({'x': x, 'y': np.interp(x, df_seabed['x'], df_seabed['y'])})
    xmin_spline, xmax_spline = spline_seabed['x'].min(), spline_seabed['x'].max()
    ymin_spline, ymax_spline = spline_seabed['y'].min(), spline_seabed['y'].max()
    node_in_lim = noeud[(noeud['x'] >= xmin_spline) & (noeud['x'] <= xmax_spline) & (noeud['y'] >= ymin_spline) & (noeud['y'] <= ymax_spline)]
    # Le dataframe node_in_lim contient les noeuds proches dans df_seabed
    node_seabed = pd.DataFrame()
    # Si la distance entre node_lin et line_seabed est dans l'intervalle de 0.1, alors node_lin est dans line_seabed
    for index, row in node_in_lim.iterrows():
        # Calculer la distance entre le noeud et la ligne
        x0, y0 = row['x'], row['y']
        dist = np.sqrt((spline_seabed['x'] - x0) ** 2 + (spline_seabed['y'] - y0) ** 2)
        if dist.min() < 0.01:
            # Ajouter au dataframe node_seabed
            node_seabed = pd.concat([node_seabed, pd.DataFrame([row])], ignore_index=True)

#    add point to freesurface
# create dataframe with two rows
    points_to_adds = pd.DataFrame({'x': [model_xmin, model_xmax,xmax-0.01], 'H': [df_freesurf['H'].iloc[0], df_freesurf['H'].iloc[-1], df_freesurf['H'].iloc[-1]]})
    # add the two rows to df_seabed
    df_freesurf = pd.concat([df_freesurf, points_to_adds], ignore_index=True)
    #interpolate df_freesurf['x','H'] to the node_seabed['x']
    node_seabed['H'] = np.interp(node_seabed['x'], df_freesurf['x'], df_freesurf['H'])
    min_H, max_H = node_seabed['H'].min(), node_seabed['H'].max()
    print('H=   min:', min_H, 'max:', max_H)
        
    size_closest_nodes = len(node_seabed)
    print('nombre de noeuds : ',size_closest_nodes)
    head_line = node_seabed[['identifiant', 'H']]
    head_line['identifiant'] = head_line['identifiant'].astype(int)
    name=bound_name+".n"
    head_line.to_csv(name,index=False,header=False,sep=' ')
    if visual:
        plt.figure(figsize=(10, 10))
        plt.plot(df_seabed['x'], df_seabed['y'], 'o',label='coord file seabed')
        #plot spline of df_seabed
        # above node in lim
 
        plt.plot(spline_seabed['x'], spline_seabed['y'], '-', label='spline_seabed')
        plt.plot(node_seabed['x'], node_seabed['y'], 'o',label='node_seabed')
        plt.plot(df_freesurf['x'], df_freesurf['H'], 'o',label='freesurface')
        plt.plot(node_seabed['x'], node_seabed['H'], 'o',label='node_seabed H')
        plt.legend()

        plt.show()  


def extract_nodes_lin_Q_vert(meshName, bound_name, x_bound, Q_top, Q_bot, z_apply):
    print(os.getcwd())
    filename = meshName + ".mail"
    try:
        file_path = find_file(filename, os.getcwd())
        print(f"Fichier {filename} trouvé dans {file_path}")
    except FileNotFoundError as e:
        print(e)
        return

    # Read the '.mail' file and extract the coordinates
    with open(file_path, 'r') as file:
        content = file.read()
        lines = content.split('\n')
        set_up = lines[0]
        nb_node_elt = set_up.split() 

    nb_node=int(nb_node_elt[0])
    nb_elt=int(nb_node_elt[1])

    # print('Domain=   nb_noeud:',nb_node_elt[0],'nb_elt:',nb_node_elt[1])
    # Read the 'job1.mail' file and extract the coordinates
    with open(file_path, 'r') as file:
        content = file.read()
        lines = content.split('\n')
        node_lines = lines[1:nb_node+1]
        coordinates = [line.split() for line in node_lines]
        x = [float(coord[0]) for coord in coordinates]
        y = [float(coord[1]) for coord in coordinates]
        noeud = pd.DataFrame({'x': x, 'y': y})

    # Reset the index and rename the column
    noeud = noeud.reset_index()
    noeud = noeud.rename(columns={'index': 'identifiant'})

    # Update the 'noeud' dataframe
    noeud['identifiant'] = noeud['identifiant'] + 1
    # find all nodes in the line x=x_bound and between y=0 and y=z_apply
# extract node
    node_in_line= noeud[(noeud['x'] == x_bound) & (noeud['y'] >=0 ) & (noeud['y'] <= z_apply)].copy()
    node_in_line['Q'] = Q_bot + (node_in_line['y'] / z_apply) * (Q_top - Q_bot)

# write in a vector node_in_line['identifiant'] node_in_line['head']
    head_line = node_in_line[['identifiant', 'Q']]
    outname=bound_name+".n"
    head_line.to_csv(outname, index=False, header=False, sep=' ')
    return node_in_line

def extract_nodes_temp(meshName,z_sed,temp_sea):
    print(os.getcwd())
    name = "test"
    filename = meshName + ".mail"
    filesave = name + ".n"
    try:
        file_path = find_file(filename, os.getcwd())
        print(f"Fichier {filename} trouvé dans {file_path}")
    except FileNotFoundError as e:
        print(e)


        # Read the '.mail' file and extract the coordinates
    with open(file_path, 'r') as file:
        content = file.read()
        lines = content.split('\n')
        set_up = lines[0]
        nb_node_elt = set_up.split() 

    nb_node=int(nb_node_elt[0])
    nb_elt=int(nb_node_elt[1])

    print('Domain=   nb_noeud:',nb_node_elt[0],'nb_elt:',nb_node_elt[1])
        # Read the 'job1.mail' file and extract the coordinates
    with open(file_path, 'r') as file:
        content = file.read()
        lines = content.split('\n')
        node_lines = lines[1:nb_node+1]
        coordinates = [line.split() for line in node_lines]
        x = [float(coord[0]) for coord in coordinates]
        y = [float(coord[1]) for coord in coordinates]
        noeud = pd.DataFrame({'x': x, 'y': y})

        # Reset the index and rename the column
    noeud = noeud.reset_index()
    noeud = noeud.rename(columns={'index': 'identifiant'})

        # Update the 'noeud' dataframe
    noeud['identifiant'] = noeud['identifiant'] + 1

    # Trouver tous les nœuds entre y=z_sed et y=z_sed+0.5
    node_in_line = noeud[(noeud['y'] >= z_sed+0.5) & (noeud['y'] <= z_sed+1.5)]

    # Supposons que vous vouliez définir une valeur pour 'T' et une autre colonne, par exemple 'head'
    node_in_line['T'] = temp_sea
    # Assurez-vous que 'head' est correctement défini ou mis à jour ici, par exemple :
    # node_in_line['head'] = valeur_appropriée

    # Écrire dans un vecteur node_in_line['identifiant'] et node_in_line['T']
    T_line = node_in_line[['identifiant', 'T']]
    T_line.to_csv(filesave, index=False, header=False, sep=' ')



# Add a single top-level definition at the end of the file

def build_gmsh_mesh_from_linestopbottom(
    x_list, z_haut_list, z_bas_list, xr1=0, xr2=0, z_zns=78.5,
    resolution_max=0.5, resolution_min=0.01,
    resolution_seabed=0.5, resolution_bound=0.5, Dist_bound=5,
    mesh_filename="profile_mesh.msh", visualize=False
):
    """
    Create a GMSH mesh from a profile defined by x, z_haut, z_bas.
    Two mesh zones: fine above z_zns, coarse below.
    Args:
        x_list: list or np.array of x coordinates
        z_haut_list: list or np.array of z_haut (top) values
        z_bas_list: list or np.array of z_bas (bottom) values
        mesh_filename: output .msh file name
        visualize: if True, show mesh in gmsh GUI

    _________________                                               __________________ z_haut
                     _______                                ________
                          |xr1______________________________         |x_2                 |L
    ____________________________________________________________________________________ 
    |                     |xr1                                       |x_2                 |L
    |                     |                                          |L                   
    |                     |xr1                                       |x_2                 |L
    _____________________________________________________________________________________ z_bot
    

    
    
    Parameters:
     ----------
    meshName : str, name of the meshd
    z_bot : float, bottom of the sediment, with a default value if not provided
    z_haut : float, elevation of the sea domain, with a default value if not provided
    resolution_max : float, maximum resolution of the mesh
    resolution_streambed : float, minimum resolution of the seabed
    resolution_bound : float, resolutions for the mesh at the boundaries in the sea domain
    Dist_bound : float, distance from the boundaries for the mesh where the mesh will become the maximum resolution
    Dist_streambed : float, distance from the seabed for the mesh where the mesh will become the maximum resolution

    
    Imposed parameters:
    ------------------
    L : float, length of the domain fixed at 2xr1+xr2
    the reference point is fixed at x=0, z=z_bot
    LcMin : This is the minimum mesh size. In the context of the "Threshold" field, it is the mesh size below the threshold. It is used to control the granularity of the mesh. A smaller value will result in a finer mesh.
    LcMax : This is the maximum mesh size. In the context of the "Threshold" field, it is the mesh size above the threshold. It is used to control the granularity of the mesh. A larger value will result in a coarser mesh.
    DistMin : the distance where the mesh size is interpolated between `LcMin` and `LcMax`.
    DistMax : the maximum distance above which the mesh size is egal to the maximum mesh size `LcMax`.
    Threshold field is used to create a transition between two prescribed mesh sizes (`LcMin` and `LcMax`) based on the distance to a given set of points, curves, surfaces or volumes (defined by the "Distance" field). The distances below and above which the mesh sizes are interpolated are defined by `DistMin` and `DistMax`.
    ----------
    Returns
    -------
    The "meshName".msh file
    
    
    """
    import gmsh
    import numpy as np
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal",0)#     
    gmsh.model.add("profile_mesh")


    # Utilise uniquement les coins extrêmes pour faire un rectangle (4 lignes, pas de surface)
    x_left = float(x_list[0])
    x_right = float(x_list[-1])
    z_top_left = float(z_haut_list[0])
    z_top_right = float(z_haut_list[-1])
    z_bot_left = float(z_bas_list[0])
    z_bot_right = float(z_bas_list[-1])
    z_top=max(z_top_left,z_top_right)
    z_bot=min(z_bot_left,z_bot_right)
    # Création des points pour le profil supérieur (z_haut_list), intermédiaire (z_zns), et inférieur (z_bas_list)
    n = len(x_list)
    pt_sup = [gmsh.model.geo.addPoint(float(x), float(z_haut), 0, resolution_min) for x, z_haut in zip(x_list, z_haut_list)]
    pt_mid = [gmsh.model.geo.addPoint(float(x), float(z_zns), 0, resolution_max) for x in x_list]
    pt_bot = [gmsh.model.geo.addPoint(float(x), float(z_bot), 0, resolution_seabed) for x in x_list]

    # Surface supérieure : entre z_haut_list et z_zns
    # Correction du sens pour le contour fermé (antihoraire)
    # Correction : cycle explicite, pas d'inversion
    spline_sup = gmsh.model.geo.addSpline(pt_sup)
    spline_mid = gmsh.model.geo.addSpline(pt_mid[::-1])  # z_zns de droite à gauche
    l_sup_left = gmsh.model.geo.addLine(pt_mid[0],pt_sup[0])  # gauche haut -> gauche z_zns
    l_sup_right = gmsh.model.geo.addLine( pt_sup[-1],pt_mid[-1])  # droite z_zns -> droite haut
    cl_sup = gmsh.model.geo.addCurveLoop([spline_sup, l_sup_right, spline_mid, l_sup_left])
    surf_sup = gmsh.model.geo.addPlaneSurface([cl_sup])

    # Surface inférieure : entre z_zns et z_bas_list
    spline_bot = gmsh.model.geo.addSpline(pt_bot[::-1])
    l_inf_left = gmsh.model.geo.addLine(pt_bot[0],pt_mid[0])  # gauche bas -> gauche z_zns
    l_inf_right = gmsh.model.geo.addLine( pt_mid[-1],pt_bot[-1])  # droite zns -> droite bas
    spline_mid2 = gmsh.model.geo.addSpline(pt_mid) # gauche  -> droite z_zns
    cl_inf = gmsh.model.geo.addCurveLoop([spline_mid2, l_inf_right, spline_bot, l_inf_left])
    surf_inf = gmsh.model.geo.addPlaneSurface([cl_inf])

    gmsh.model.geo.synchronize()

    # Définir la taille du maillage pour chaque ligne de points
    gmsh.model.mesh.setSize([(0, tag) for tag in pt_sup], resolution_min)
    gmsh.model.mesh.setSize([(0, tag) for tag in pt_mid], resolution_max)
    gmsh.model.mesh.setSize([(0, tag) for tag in pt_bot], resolution_seabed)

    # Générer le maillage
    gmsh.model.mesh.generate(2)

    # Sauvegarder le maillage
    gmsh.write(mesh_filename)

    if visualize:
        gmsh.fltk.run()
    gmsh.finalize()


import meshio
import numpy as np

def remove_elements_above_z_haut(mesh_filename, x_list, z_haut_list, output_filename=None):
    """
    Supprime les éléments dont tous les sommets sont au-dessus de la ligne z_haut(x).
    """
    import meshio
    import numpy as np

    mesh = meshio.read(mesh_filename)
    points = mesh.points
    cells = mesh.cells_dict["triangle"]

    # Interpolation de z_haut(x) pour chaque x
    z_haut_interp = lambda x: np.interp(x, x_list, z_haut_list)

    def is_below_z_haut(cell):
        pts = points[cell]
        return np.all(pts[:, 1] <= z_haut_interp(pts[:, 0]))

    mask = np.array([is_below_z_haut(cell) for cell in cells])
    filtered_cells = cells[mask]

    # Correction ici : ne garder que les cell_data des triangles
    filtered_cell_data = {}
    if mesh.cell_data:
        for key, value in mesh.cell_data.items():
            # value est une liste, chaque entrée correspond à un type de cellule
            # On cherche l'entrée correspondant aux triangles
            filtered_value = []
            for cell_block, cell_array in zip(mesh.cells, value):
                if cell_block.type == "triangle":
                    filtered_value.append(np.array(cell_array)[mask])
                else:
                    filtered_value.append(cell_array)
            filtered_cell_data[key] = filtered_value

    mesh_filtered = meshio.Mesh(
        points=points,
        cells=[("triangle", filtered_cells)],
        point_data=mesh.point_data,
        cell_data=filtered_cell_data if filtered_cell_data else None
    )

    if output_filename is None:
        output_filename = mesh_filename.replace(".msh", "_filtered.msh")
    meshio.write(output_filename, mesh_filtered)
    print(f"Maillage filtré sauvegardé sous {output_filename}")



def extract_nodes_lin_H_vert(meshName, bound_name, x_bound, head_top, head_bot, z_apply):
    # print(os.getcwd())
    filename = meshName + ".mail"

    noeud = read_mail_file(filename)

    # Update the 'noeud' dataframe
    noeud['identifiant'] = noeud['identifiant'] + 1
    # find all nodes in the line x=x_bound and between y=0 and y=z_apply
# extract node
    node_in_line= noeud[(noeud['x'] == x_bound) & (noeud['y'] >= 0) & (noeud['y'] <= z_apply)].copy()
    node_in_line['head'] = head_bot + (node_in_line['y'] / z_apply) * (head_top - head_bot)

# write in a vector node_in_line['identifiant'] node_in_line['head']
    head_line = node_in_line[['identifiant', 'head']]
    outname=bound_name+".n"
    head_line.to_csv(outname, index=False, header=False, sep=' ')
    return node_in_line

    
def extract_nodes_lin_H_seabed(meshName,name,z_sed, head_left=12.8, head_right=17.17,L=200):
    # print(os.getcwd())
    filename = meshName + ".mail"
    filesave=name+".n"
    # Update the 'noeud' dataframe
    noeud=read_mail_file(meshName)
    y_line_seased= z_sed
    # find all nodes in the line y=0
    node_in_line = noeud[noeud['y'] == y_line_seased].copy()   
    # calculate head of each node of line_seased
    #head(0)=14 head(200)=15 head(x)=head(0)+(head(200)-head(0))/(200-0)*x
    node_in_line['head'] = head_left + (head_right - head_left) / (L - 0) * node_in_line['x']

    # write in a vector node_in_line['identifiant'] node_in_line['head']

    head_line = node_in_line[['identifiant', 'head']]
    head_line.to_csv(filesave, index=False, header=False, sep=' ')
    return node_in_line


def extract_nodes_cable(meshName,name,cable_x,cable_y,cable_r):
    # print(os.getcwd())
    filename = meshName + ".mail"
    filesave = name+".n"
    noeud=read_mail_file(meshName)

    # Ajoutez une nouvelle colonne 'distance_to_center' qui calcule la distance de chaque nœud au centre du cercle
    noeud['distance_to_center'] = np.sqrt((noeud['x'] - cable_x)**2 + (noeud['y'] - cable_y)**2)
    #find in noeud all nodes in the squate cable_x-cable_r<=x<=cable_x+cable_r and cable_y-cable_r<=y<=cable_y+cable_r
    noeud_sq=noeud[(noeud['x'] >= cable_x - cable_r-1) & (noeud['x'] <= cable_x + cable_r+1) & (noeud['y'] >= cable_y - cable_r-1 )& (noeud['y'] <= cable_y + cable_r+1)].copy()
    #order noeud by distance to center
    noeud = noeud.sort_values(by='distance_to_center')

    # Filtrez les nœuds dont la distance au centre est égale au rayon du cercle

    closest_nodes = noeud[np.isclose(noeud['distance_to_center'], cable_r, atol=1e-1)]



    size_closest_nodes = len(closest_nodes)
    # print('nombre de noeuds sur le cable: ',size_closest_nodes)
    nodes_cable=pd.DataFrame({"id":closest_nodes["identifiant"]})
    nodes_cable["no"]="no"
    nodes_cable=nodes_cable[nodes_cable.columns[::-1]]

    nodes_cable.to_csv(filesave,index=False,header=False,sep=' ')
    # print('fichier cable.n créé')

def read_dat2mesh(meshName,datFile,z_bot=-100,z_sea=15,cable_y=-3, cable_r=0.246/2,
                  resolution_max = 10,resolution_cable = 0.01,
                  resolution_seabed=0.5,resolution_bound=0.5,Dist_bound=5,Dist_seabed=5,Dist_cable=50,visual=True):


    '''
    Creates a flat mesh for a simulation. The discretisation is thinner around the cable, the seabed and the boundaries of the sea domain.
    
    Figure of parameters:
     ----------
    _____________________________________________________________________________________ z_sea
    
    |                     |xr1                                       |x_2                 |L
    ____________________________________________________________________________________ z_sed
    |                     |xr1                                       |x_2                 |L
    |                     |                     o(cable_x,cable_y,cable_r)                |L                   
    |                     |xr1                                       |x_2                 |L
    _____________________________________________________________________________________ z_bot
    
    Figure of Discretisation:
     ----------
    _____________________________________________________________________________________ z_sea
    
    | resolution_bound    |x_1                                       |x_2 resolution_bound|L
    __________________resolution_seabed___________________________________________________ z_sed
    |                     |x_1                                       |x_2                 |L
    |                     |                     o(resolution_cable)  |x_2                 |L                   
    |                     |x_1                                       |x_2                 |L
    _____________________________________________________________________________________ z_bot
    
    
    Parameters:
     ----------
    meshName : str, name of the meshd
    z_bot : float, bottom of the sediment, with a default value if not provided
    z_sea : float, elevation of the sea domain, with a default value if not provided
    cable_y : float, coordinates of the cable, with default values if not provided
    cable_r : float, radius of the cable, with a default value if not provided
    resolution_max : float, maximum resolution of the mesh
    resolution_cable : float, minimun resolution of the cable
    resolution_seabed : float, minimum resolution of the seabed
    resolution_bound : float, resolutions for the mesh at the boundaries in the sea domain
    Dist_bound : float, distance from the boundaries for the mesh where the mesh will become the maximum resolution
    Dist_seabed : float, distance from the seabed for the mesh where the mesh will become the maximum resolution
    Dist_cable : float, distances for the mesh
    
    Imposed parameters:
    ------------------
    cable_x : float, x-coordinate of the cable fixed at the middle of the domain
    L : float, length of the domain fixed at 2xr1+xr2
    the reference point is fixed at x=0, z=z_bot
    LcMin : This is the minimum mesh size. In the context of the "Threshold" field, it is the mesh size below the threshold. It is used to control the granularity of the mesh. A smaller value will result in a finer mesh.
    LcMax : This is the maximum mesh size. In the context of the "Threshold" field, it is the mesh size above the threshold. It is used to control the granularity of the mesh. A larger value will result in a coarser mesh.
    DistMin : the distance where the mesh size is interpolated between `LcMin` and `LcMax`.
    DistMax : the maximum distance above which the mesh size is egal to the maximum mesh size `LcMax`.
    Threshold field is used to create a transition between two prescribed mesh sizes (`LcMin` and `LcMax`) based on the distance to a given set of points, curves, surfaces or volumes (defined by the "Distance" field). The distances below and above which the mesh sizes are interpolated are defined by `DistMin` and `DistMax`.
    ----------
    Returns
    -------
    The "meshName".msh file
    
    '''

    import gmsh
    import sys
    import math


    # Initialisation de Gmsh
    gmsh.initialize()
    model = gmsh.model
    factory = model.occ




    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", resolution_cable)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", resolution_max)


# read the data file and plot the line
    with open(datFile, 'r') as f:
        lines = f.readlines()
        x = []
        y = []
        for line in lines:
            line = line.split()
            x.append(float(line[0]))
            y.append(float(line[1]))
    # create dataframe
    df = pd.DataFrame({'x': x, 'y':y})
    # order by x
    df = df.sort_values(by='x')
    x_min = min(x)
    x_max = max(x)
    y_min = min(y)
    y_max = max(y)

    xr1 = x_min
    xr2 = x_max

    z_sed1 = df['y'].values[0]
    z_sed2 = df['y'].values[-1]
    h_sed1 = abs(z_sed1 - z_bot)
    h_sed2 = abs(z_sed2 - z_bot)
    L = xr1 + xr2
    cable_x = L / 2
    lc = resolution_max

    points_ids = []
    with open(datFile) as fin:
        for line in fin:
            try:
                x, y = map(float, line.split())
                pid = factory.addPoint(x, y, 0)
                points_ids.append(pid)
            except ValueError:
                print("Skipped line: ", line)

    if len(points_ids) > 1:
        spline = factory.addBSpline(points_ids)
        # print("Spline créée.")
    else:
        print("Pas assez de points pour créer une spline.")
        spline = None


    factory.synchronize()



    p1=points_ids[0]
    p2=points_ids[-1]
    p3=factory.addPoint(xr2, z_bot, 0)
    p4=factory.addPoint(xr1, z_bot, 0)
    factory.synchronize()
    l2=factory.addLine(p3, p2)
    l3=factory.addLine(p3, p4)
    l4=factory.addLine(p4, p1)
    factory.synchronize()
    Cloop1=factory.addCurveLoop([l4, spline, -l2, l3])
    factory.synchronize()
    Psurf1=factory.addPlaneSurface([Cloop1])
    factory.synchronize()
    PhyGp1=model.addPhysicalGroup(2, [Psurf1])
    factory.synchronize()



        

    
    tag_point = len(model.getEntities(0))
    tag_line = len(model.getEntities(1))
    tag_surface = len(model.getEntities(2))
    tag_physical = len(model.getEntities(3))
  #  print(tag_point, tag_line, tag_surface, tag_physical)
    p5=factory.addPoint(xr2, z_sea, 0)
    factory.synchronize()    
    p6=factory.addPoint(xr1, z_sea, 0)
    factory.synchronize()
    l5=factory.addLine(p5, p2)
    factory.synchronize()    
    l6=factory.addLine(p5, p6)
    l7=factory.addLine(p6, p1)
    factory.synchronize()
    Cloop2=factory.addCurveLoop([l7, spline, -l5, l6])
    factory.synchronize()
    Psurf2=factory.addPlaneSurface([Cloop2])
    factory.synchronize()
    PhyGp2=model.addPhysicalGroup(2, [Psurf2])
    factory.synchronize()       

    p7 = factory.addPoint(0, z_sed1, 0)
    p8 = factory.addPoint(0, z_bot, 0)
    factory.synchronize()    
    l8 = factory.addLine(p7, p1)
    l9 = factory.addLine(p4, p8)
    l10 = factory.addLine(p8, p7)
    factory.synchronize()    
    Cloop3=factory.addCurveLoop([l8, -l4, l9, l10])
    factory.synchronize()     
    Psurf3=factory.addPlaneSurface([Cloop3])
    factory.synchronize()    
    PhyGp3=model.addPhysicalGroup(2, [Psurf3])
    factory.synchronize()

    p9 = factory.addPoint(0, z_sea, 0)
    factory.synchronize()
    l11 = factory.addLine(p9, p6)
    l12 = factory.addLine(p7, p9)
    factory.synchronize()    
    Cloop4=factory.addCurveLoop([l11, l7, -l8, l12])
    factory.synchronize()    
    Psurf4=factory.addPlaneSurface([Cloop4])
    factory.synchronize()    
    PhyGp4=model.addPhysicalGroup(2, [Psurf4])
    factory.synchronize()
        
    p10 = factory.addPoint(L, z_sed2, 0)
    p12 = factory.addPoint(L, z_bot, 0)
    factory.synchronize()    
    l13 = factory.addLine(p10, p2)
    l14 = factory.addLine(p3, p12)
    l15 = factory.addLine(p12, p10)
    factory.synchronize()    
    Cloop5=factory.addCurveLoop([l13, -l2, l14, l15])
    factory.synchronize()
    Psurf5=factory.addPlaneSurface([Cloop5])
    factory.synchronize()    
    PhyGp5=model.addPhysicalGroup(2, [Psurf5])
    factory.synchronize()

        
    p13 = factory.addPoint(L, z_sea, 0)
    factory.synchronize()
    l16 = factory.addLine(p13, p5)
    l17 = factory.addLine(p10, p13)
    factory.synchronize()
    Cloop6=factory.addCurveLoop([l16, l5, -l13, l17])
    factory.synchronize()    
    Psurf6=factory.addPlaneSurface([Cloop6])
    factory.synchronize()    
    PhyGp6=model.addPhysicalGroup(2, [Psurf6])    
    factory.synchronize()
    # Ajouter le cercle au modèle
    x_circle = cable_x  # Utilisez la variable cable_x pour la position x du cercle
    y_circle = cable_y  # Utilisez la variable cable_y pour la position y du cercle
    z_circle = 0  # La position z du cercle, ajustez selon le besoin
    r_circle = cable_r  # Le rayon du cercle, ici utilisé cable_r
    pcircle = factory.addPoint(cable_x, cable_y, z_circle, resolution_cable)
    factory.synchronize()
    circle = factory.addDisk(x_circle, y_circle, z_circle, r_circle, r_circle)
    factory.synchronize()
    # Fragmenter le rectangle et le cercle pour créer un trou
    factory.cut([(2,Psurf1)], [(2, circle)])
    factory.synchronize()
    # update the physical group
    model.addPhysicalGroup(2, [Psurf1])
 # Fragmenter les entités
    factory.fragment([(2, Cloop1), (2, Cloop2),(2,Cloop3),(2,Cloop4),(2,Cloop5),(2,Cloop6)], [(2, Cloop1), (2, Cloop2),(2,Cloop3),(2,Cloop4),(2,Cloop5),(2,Cloop6)])
    factory.synchronize()

    # Ajouter le premier point au modèle
    point1 = factory.addPoint(cable_x, 1, 0)

    # Ajouter le deuxième point au modèle
    point2 = factory.addPoint(cable_x, cable_y-2*cable_r, 0)

    # Ajouter la ligne au modèle
    line = factory.addLine(point1, point2)

    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter un champ "Distance"
    field_id = model.mesh.field.add("Distance")

    # Distance ligne
    model.mesh.field.setNumbers(field_id, "CurvesList", [line])

    # Ajouter un champ "Threshold"
    threshold_cable = model.mesh.field.add("Threshold")

    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_cable, "IField", field_id)
    model.mesh.field.setNumber(threshold_cable, "LcMin", resolution_cable)
    model.mesh.field.setNumber(threshold_cable, "LcMax",resolution_max)
    model.mesh.field.setNumber(threshold_cable, "DistMin", cable_r+0.01)
    model.mesh.field.setNumber(threshold_cable, "DistMax", Dist_cable)  # Distance entre x=50 et le centre

     # meshing
    # Ajouter un champ "Distance"
    field_id1 = model.mesh.field.add("Distance")
    # Distance ligne gauche
    model.mesh.field.setNumbers(field_id1, "EdgesList", [l12])
    factory.synchronize()
    # Ajouter un champ "Threshold"
    threshold_id1 = model.mesh.field.add("Threshold")
    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_id1, "IField", field_id1)
    model.mesh.field.setNumber(threshold_id1, "LcMin", resolution_bound)
    model.mesh.field.setNumber(threshold_id1, "LcMax", resolution_max)
    model.mesh.field.setNumber(threshold_id1, "DistMin", 0)
    model.mesh.field.setNumber(threshold_id1, "DistMax", Dist_bound)
    factory.synchronize()

    # Ajouter un champ "Distance"
    field_id2 = model.mesh.field.add("Distance")
    # Distance ligne droite
    model.mesh.field.setNumbers(field_id2, "EdgesList", [l17])
    # Ajouter un champ "Threshold"
    threshold_id2 = model.mesh.field.add("Threshold")
    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_id2, "IField", field_id2)
    model.mesh.field.setNumber(threshold_id2, "LcMin", resolution_bound)
    model.mesh.field.setNumber(threshold_id2, "LcMax", resolution_max)
    model.mesh.field.setNumber(threshold_id2, "DistMin", 0)
    model.mesh.field.setNumber(threshold_id2, "DistMax", Dist_bound)
    factory.synchronize()
 

    # Ajouter un champ "Distance"
    field_id3 = model.mesh.field.add("Distance")
    # Distance ligne
    model.mesh.field.setNumbers(field_id3, "EdgesList", [spline])
    # Ajouter un champ "Threshold"
    threshold_id3 = model.mesh.field.add("Threshold")
    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_id3, "IField", field_id3)
    model.mesh.field.setNumber(threshold_id3, "LcMin", resolution_seabed)
    model.mesh.field.setNumber(threshold_id3, "LcMax", resolution_max)
    model.mesh.field.setNumber(threshold_id3, "DistMin", 0)
    model.mesh.field.setNumber(threshold_id3, "DistMax", Dist_seabed)
    factory.synchronize()   

    nb_nodes_seabed_l=int(xr1/resolution_seabed)
    # print('nb_nodes_seabed_l:',nb_nodes_seabed_l)    
      # Ajouter un champ "Distance"
    field_id4 = model.mesh.field.add("Distance")
    # Distance ligne droite
    model.mesh.field.setNumbers(field_id4, "EdgesList", [l8])
    model.mesh.field.setNumbers(field_id4, "PointsList", [p1])
    # Ajouter un champ "Threshold"
    threshold_id4 = model.mesh.field.add("Threshold")
    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_id4, "IField", field_id4)
    model.mesh.field.setNumber(threshold_id4, "LcMin", resolution_seabed)
    model.mesh.field.setNumber(threshold_id4, "LcMax", resolution_max)
    model.mesh.field.setNumber(threshold_id4, "DistMin", 0)
    model.mesh.field.setNumber(threshold_id4, "DistMax", Dist_seabed)
    factory.synchronize()
    
    nb_nodes_seabed_r=int((L-xr2)/resolution_seabed)
    # print('nb_nodes_seabed_r:',nb_nodes_seabed_r)    
      # Ajouter un champ "Distance"
    field_id5 = model.mesh.field.add("Distance")
    # Distance ligne droite
    model.mesh.field.setNumbers(field_id5, "EdgesList", [l13])
    model.mesh.field.setNumbers(field_id5, "PointsList", [p2])
    # Ajouter un champ "Threshold"
    threshold_id5 = model.mesh.field.add("Threshold")
    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_id5, "IField", field_id5)
    model.mesh.field.setNumber(threshold_id5, "LcMin", resolution_seabed)
    model.mesh.field.setNumber(threshold_id5, "LcMax", resolution_max)
    model.mesh.field.setNumber(threshold_id5, "DistMin", 0)
    model.mesh.field.setNumber(threshold_id5, "DistMax", Dist_seabed)
    factory.synchronize()  


            # Ajouter un champ "Min"
    min_field_id = model.mesh.field.add("Min")

        # Définir les paramètres pour le champ "Min"
    model.mesh.field.setNumbers(min_field_id, "FieldsList", [threshold_cable,threshold_id1,threshold_id2,threshold_id3,threshold_id4,threshold_id5])
        # Définir le champ "Min" comme champ de taille de maillage
    model.mesh.field.setAsBackgroundMesh(min_field_id)

        # Synchroniser le modèle
    factory.synchronize()



    # Génération du maillage
    model.mesh.generate(2)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    factory.synchronize()
    # Mesh save
    meshName_msh = meshName + ".msh"
    # print('nom fichier:',meshName_msh)
    gmsh.write(meshName_msh)
    if visual:
        gmsh.fltk.run()
    gmsh.finalize()
    return cable_x,L   
        
        
def meshplat(meshName,L=200, xr1=50,xr2=150,z_bot=-100, z_sea=15, z_sed=0, cable_x=100, cable_y=-3, 
             cable_r=0.246/2,resolution_max = 10, resolution_cable = 0.01,resolution_seabed=0.5,
             resolution_bound=0.5,Dist_bound=5, Dist_seabed=5,Dist_cable=50,visual=True):
    '''
    Creates a flat mesh for a simulation. The discretisation is thinner around the cable, the seabed and the boundaries of the sea domain.
    
    Figure of parameters:
     ----------
    _____________________________________________________________________________________ z_sea
    
    |                     |xr1                                       |x_2                 |L
    ____________________________________________________________________________________ z_sed
    |                     |xr1                                       |x_2                 |L
    |                     |                     o(cable_x,cable_y,cable_r)                |L                   
    |                     |xr1                                       |x_2                 |L
    _____________________________________________________________________________________ z_bot
    
    Figure of Discretisation:
     ----------
    _____________________________________________________________________________________ z_sea
    
    | resolution_bound    |x_1                                       |x_2 resolution_bound|L
    __________________resolution_seabed___________________________________________________ z_sed
    |                     |x_1                                       |x_2                 |L
    |                     |                     o(resolution_cable)  |x_2                 |L                   
    |                     |x_1                                       |x_2                 |L
    _____________________________________________________________________________________ z_bot
    
    
    Parameters:
     ----------
    meshName : str, name of the mesh
    L : float, length of the domain for the simulation, with a default value if not provided
    xr1 : float, coordinate first rectangle, with a default value if not provided
    xr2 : float, coordinate last rectangle, with a default value if not provided
    z_bot : float, bottom of the sediment, with a default value if not provided
    z_sea : float, elevation of the sea domain, with a default value if not provided
    z_sed : float, elevation of the seabed, with a default value if not provided
    cable_x, cable_y : float, coordinates of the cable, with default values if not provided
    cable_r : float, radius of the cable, with a default value if not provided
    LcMin : This is the minimum mesh size. In the context of the "Threshold" field, it is the mesh size below the threshold. It is used to control the granularity of the mesh. A smaller value will result in a finer mesh.
    LcMax : This is the maximum mesh size. In the context of the "Threshold" field, it is the mesh size above the threshold. It is used to control the granularity of the mesh. A larger value will result in a coarser mesh.
    DistMin : the distance where the mesh size is interpolated between `LcMin` and `LcMax`.
    DistMax : the maximum distance above which the mesh size is egal to the maximum mesh size `LcMax`.
    Threshold field is used to create a transition between two prescribed mesh sizes (`LcMin` and `LcMax`) based on the distance to a given set of points, curves, surfaces or volumes (defined by the "Distance" field). The distances below and above which the mesh sizes are interpolated are defined by `DistMin` and `DistMax`.
    resolution_max : float, maximum resolution of the mesh
    resolution_cable : float, minimun resolution of the cable
    resolution_seabed : float, minimum resolution of the seabed
    resolution_bound : float, resolutions for the mesh at the boundaries in the sea domain
    Dist_bound : float, distance from the boundaries for the mesh where the mesh will become the maximum resolution
    Dist_seabed : float, distance from the seabed for the mesh where the mesh will become the maximum resolution
    Dist_cable : float, distances for the mesh
    ----------
    Returns
    -------
    The "meshName".msh file
    
    '''
    # Initialisation de Gmsh
    gmsh.initialize()
    
    gmsh.option.setNumber("General.Terminal",0)# 
    
    model = gmsh.model
    factory = model.occ
    mesh = model.mesh

    resolution_min = resolution_cable
    

    lc=resolution_max


    # code
    h_sed=abs(z_sed-z_bot)

    gmsh.initialize()
    model.add("test")

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", resolution_cable)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", resolution_max)

    # utilisée pour ajouter un rectangle à votre modèle dans Gmsh. Les arguments de cette fonction sont les coordonnées du coin inférieur gauche du rectangle et les longueurs des côtés du rectangle.
    # Créer les rectangles
    rect1 = factory.addRectangle(0, z_bot, 0, xr1, h_sed)
    rect3 = factory.addRectangle(xr1, z_bot, 0, xr2, h_sed)
    rect5 = factory.addRectangle(xr2, z_bot, 0, L-xr2, h_sed)

    # Synchroniser le modèle
    factory.synchronize()

    # Coordonnées du centre du cercle
    x, y, z = cable_x, cable_y, 0
    # Rayon du cercle
    r = cable_r
    # Ajouter le centre du cercle au modèle
    center = factory.addPoint(x, y, z)

    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter le cercle au modèle
    circle = factory.addDisk(x, y, z, r, r)

    # Fragmenter le rectangle et le cercle pour créer un trou
    factory.cut([(2, rect3)], [(2, circle)])
    factory.synchronize()

    # Fragmenter les entités
    factory.fragment([(2, rect1), (2, rect3),  (2, rect5)], [(2, rect1),  (2, rect3), (2, rect5)])

    # Synchroniser le modèle
    factory.synchronize()
    # Ajouter le premier point au modèle
    point1 = factory.addPoint(0, 0, 0)

    # Ajouter le deuxième point au modèle
    point2 = factory.addPoint(0, z_sea, 0)

    # Ajouter la ligne au modèle
    line = factory.addLine(point1, point2)

    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter un champ "Distance"
    field_id1 = model.mesh.field.add("Distance")

    # Distance ligne
    model.mesh.field.setNumbers(field_id1, "EdgesList", [line])
    model.mesh.field.setNumbers(field_id1, "PointsList", [point1])
    model.mesh.field.setNumbers(field_id1, "PointsList", [point2])

    # Ajouter un champ "Threshold"
    threshold_id1 = model.mesh.field.add("Threshold")

    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_id1, "IField", field_id1)
    model.mesh.field.setNumber(threshold_id1, "LcMin", resolution_bound)
    model.mesh.field.setNumber(threshold_id1, "LcMax", resolution_max)
    model.mesh.field.setNumber(threshold_id1, "DistMin", 0)
    model.mesh.field.setNumber(threshold_id1, "DistMax", Dist_bound)
    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter le premier point au modèle
    point3 = factory.addPoint(L, 0, 0)

    # Ajouter le deuxième point au modèle
    point4 = factory.addPoint(L, z_sea, 0)

    # Ajouter la ligne au modèle
    line2 = factory.addLine(point3, point4)

    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter un champ "Distance"
    field_id2 = model.mesh.field.add("Distance")

    # Distance ligne
    model.mesh.field.setNumbers(field_id2, "EdgesList", [line2])
    model.mesh.field.setNumbers(field_id2, "PointsList", [point3])
    model.mesh.field.setNumbers(field_id2, "PointsList", [point4])

    # Ajouter un champ "Threshold"
    threshold_id2 = model.mesh.field.add("Threshold")

    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_id2, "IField", field_id2)
    model.mesh.field.setNumber(threshold_id2, "LcMin", resolution_bound)
    model.mesh.field.setNumber(threshold_id2, "LcMax", resolution_max)
    model.mesh.field.setNumber(threshold_id2, "DistMin", 0)
    model.mesh.field.setNumber(threshold_id2, "DistMax", Dist_bound)
    # Synchroniser le modèle
    factory.synchronize()





    # Ajouter le premier point au modèle
    point5 = factory.addPoint(cable_x, 0, 0)

    # Ajouter le deuxième point au modèle
    point6 = factory.addPoint(cable_x, cable_y-2*cable_r, 0)

    # Ajouter la ligne au modèle
    line = factory.addLine(point5, point6)

    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter un champ "Distance"
    field_id = model.mesh.field.add("Distance")

    # Distance ligne
    model.mesh.field.setNumbers(field_id, "CurvesList", [line])

    # Ajouter un champ "Threshold"
    threshold_cable = model.mesh.field.add("Threshold")

    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_cable, "IField", field_id)
    model.mesh.field.setNumber(threshold_cable, "LcMin", resolution_cable)
    model.mesh.field.setNumber(threshold_cable, "LcMax",resolution_max)
    model.mesh.field.setNumber(threshold_cable, "DistMin", cable_r+0.01)
    model.mesh.field.setNumber(threshold_cable, "DistMax", Dist_cable)  # Distance entre x=50 et le centre
    factory.synchronize()

 # create slipe each 1 m between point1 (0,0) xr1, xr2 and point2 (L,0)
    points_ids = []   
    # Calculer et ajouter des points intermédiaires
    for x in range(0, int(L+1)):
        pid=factory.addPoint(x, z_sed, 0)
        points_ids.append(pid)
    factory.synchronize()

    if len(points_ids) > 1:
        spline = factory.addBSpline(points_ids)

    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter un champ "Distance"
# Ajouter un champ "Distance"
    field_id4 = model.mesh.field.add("Distance")
# Distance ligne
    model.mesh.field.setNumbers(field_id4, "PointsList", points_ids)  # Correction ici
    # Ajouter un champ "Threshold"
    threshold_id4 = model.mesh.field.add("Threshold")
    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_id4, "IField", field_id4)
    model.mesh.field.setNumber(threshold_id4, "LcMin", resolution_seabed)
    model.mesh.field.setNumber(threshold_id4, "LcMax", resolution_max)
    model.mesh.field.setNumber(threshold_id4, "DistMin", 0)
    model.mesh.field.setNumber(threshold_id4, "DistMax", Dist_seabed)
    factory.synchronize()
    
    
             # Ajouter un champ "Min"
    min_field_id = model.mesh.field.add("Min")

        # Définir les paramètres pour le champ "Min"
    model.mesh.field.setNumbers(min_field_id, "FieldsList", [threshold_cable,threshold_id4])
        # Définir le champ "Min" comme champ de taille de maillage
    model.mesh.field.setAsBackgroundMesh(min_field_id)

        # Synchroniser le modèle
    factory.synchronize()



    # Générer le maillage
    model.mesh.generate(2)





    # MshFileVersion settings : do not change !
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        
    # Mesh save
    meshName_msh = meshName + ".msh"
    gmsh.write(meshName_msh)
    #visua
    if visual:
        gmsh.fltk.run()
    gmsh.finalize()
    return





def meshsed(meshName,L=200, xr1=50,xr2=150,z_bot=-100, z_sed=0, cable_x=100, cable_y=-3, 
             cable_r=0.246/2,resolution_max = 10, resolution_cable = 0.01,resolution_seabed=0.5,
             resolution_bound=0.5,Dist_bound=5, Dist_seabed=5,Dist_cable=50,visual=True):
    '''
    Creates a flat mesh for a simulation. The discretisation is thinner around the cable, the seabed and the boundaries of the sea domain.
    
    Figure of parameters:
     ----------
   
    ____________________________________________________________________________________ z_sed
    |                     |xr1                                       |x_2                 |L
    |                     |                     o(cable_x,cable_y,cable_r)                |L                   
    |                     |xr1                                       |x_2                 |L
    _____________________________________________________________________________________ z_bot
    
    Figure of Discretisation:
     ----------
    _____________________________________________________________________________________ z_sea
    
    | resolution_bound    |x_1                                       |x_2 resolution_bound|L
    __________________resolution_seabed___________________________________________________ z_sed
    |                     |x_1                                       |x_2                 |L
    |                     |                     o(resolution_cable)  |x_2                 |L                   
    |                     |x_1                                       |x_2                 |L
    _____________________________________________________________________________________ z_bot

    
    Parameters:
     ----------
    meshName : str, name of the mesh
    L : float, length of the domain for the simulation, with a default value if not provided
    xr1 : float, coordinate first rectangle, with a default value if not provided
    xr2 : float, coordinate last rectangle, with a default value if not provided
    z_bot : float, bottom of the sediment, with a default value if not provided
    z_sed : float, elevation of the seabed, with a default value if not provided
    cable_x, cable_y : float, coordinates of the cable, with default values if not provided
    cable_r : float, radius of the cable, with a default value if not provided
    LcMin : This is the minimum mesh size. In the context of the "Threshold" field, it is the mesh size below the threshold. It is used to control the granularity of the mesh. A smaller value will result in a finer mesh.
    LcMax : This is the maximum mesh size. In the context of the "Threshold" field, it is the mesh size above the threshold. It is used to control the granularity of the mesh. A larger value will result in a coarser mesh.
    DistMin : the distance where the mesh size is interpolated between `LcMin` and `LcMax`.
    DistMax : the maximum distance above which the mesh size is egal to the maximum mesh size `LcMax`.
    Threshold field is used to create a transition between two prescribed mesh sizes (`LcMin` and `LcMax`) based on the distance to a given set of points, curves, surfaces or volumes (defined by the "Distance" field). The distances below and above which the mesh sizes are interpolated are defined by `DistMin` and `DistMax`.
    resolution_max : float, maximum resolution of the mesh
    resolution_cable : float, minimun resolution of the cable
    resolution_seabed : float, minimum resolution of the seabed
    resolution_bound : float, resolutions for the mesh at the boundaries in the sea domain
    Dist_bound : float, distance from the boundaries for the mesh where the mesh will become the maximum resolution
    Dist_seabed : float, distance from the seabed for the mesh where the mesh will become the maximum resolution
    Dist_cable : float, distances for the mesh
    ----------
    Returns
    -------
    The "meshName".msh file
    
    '''
    # Initialisation de Gmsh
    gmsh.initialize()
    
    gmsh.option.setNumber("General.Terminal",0)# 
    
    model = gmsh.model
    factory = model.occ
    mesh = model.mesh

    resolution_min = resolution_cable
    

    lc=resolution_max


    # code
    h_sed=abs(z_sed-z_bot)

    gmsh.initialize()
    model.add("test")

    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", resolution_cable)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", resolution_max)

    # utilisée pour ajouter un rectangle à votre modèle dans Gmsh. Les arguments de cette fonction sont les coordonnées du coin inférieur gauche du rectangle et les longueurs des côtés du rectangle.
    # Créer les rectangles
    rect1 = factory.addRectangle(0, z_bot, 0, xr1, h_sed)
    rect3 = factory.addRectangle(xr1, z_bot, 0, xr2, h_sed)
    rect5 = factory.addRectangle(xr2, z_bot, 0, L-xr2, h_sed)

    # Synchroniser le modèle
    factory.synchronize()

    # Coordonnées du centre du cercle
    x, y, z = cable_x, cable_y, 0
    # Rayon du cercle
    r = cable_r
    # Ajouter le centre du cercle au modèle
    center = factory.addPoint(x, y, z)

    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter le cercle au modèle
    circle = factory.addDisk(x, y, z, r, r)

    # Fragmenter le rectangle et le cercle pour créer un trou
    factory.cut([(2, rect3)], [(2, circle)])
    factory.synchronize()

    # Fragmenter les entités
    factory.fragment([(2, rect1),(2, rect3),  (2, rect5)], [(2, rect1),  (2, rect3), (2, rect5)])

    # Synchroniser le modèle
    factory.synchronize()
    # Ajouter le premier point au modèle
    point1 = factory.addPoint(0, 0, 0)

    # Ajouter le premier point au modèle
    point3 = factory.addPoint(L, 0, 0)

    factory.synchronize()





    # Ajouter le premier point au modèle
    point5 = factory.addPoint(cable_x, 0, 0)

    # Ajouter le deuxième point au modèle
    point6 = factory.addPoint(cable_x, cable_y-2*cable_r, 0)

    # Ajouter la ligne au modèle
    line = factory.addLine(point5, point6)

    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter un champ "Distance"
    field_id = model.mesh.field.add("Distance")

    # Distance ligne
    model.mesh.field.setNumbers(field_id, "CurvesList", [line])

    # Ajouter un champ "Threshold"
    threshold_cable = model.mesh.field.add("Threshold")

    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_cable, "IField", field_id)
    model.mesh.field.setNumber(threshold_cable, "LcMin", resolution_cable)
    model.mesh.field.setNumber(threshold_cable, "LcMax",resolution_max)
    model.mesh.field.setNumber(threshold_cable, "DistMin", cable_r+0.01)
    model.mesh.field.setNumber(threshold_cable, "DistMax", Dist_cable)  # Distance entre x=50 et le centre
    factory.synchronize()

 # create slipe each 1 m between point1 (0,0) xr1, xr2 and point2 (L,0)
    points_ids = []   
    # Calculer et ajouter des points intermédiaires
    for x in range(0, int(L+1)):
        pid=factory.addPoint(x, z_sed, 0)
        points_ids.append(pid)
    factory.synchronize()

    if len(points_ids) > 1:
        spline = factory.addBSpline(points_ids)

    # Synchroniser le modèle
    factory.synchronize()

    # Ajouter un champ "Distance"
# Ajouter un champ "Distance"
    field_id4 = model.mesh.field.add("Distance")
# Distance ligne
    model.mesh.field.setNumbers(field_id4, "PointsList", points_ids)  # Correction ici
    # Ajouter un champ "Threshold"
    threshold_id4 = model.mesh.field.add("Threshold")
    # Définir les paramètres pour le champ "Threshold"
    model.mesh.field.setNumber(threshold_id4, "IField", field_id4)
    model.mesh.field.setNumber(threshold_id4, "LcMin", resolution_seabed)
    model.mesh.field.setNumber(threshold_id4, "LcMax", resolution_max)
    model.mesh.field.setNumber(threshold_id4, "DistMin", 0)
    model.mesh.field.setNumber(threshold_id4, "DistMax", Dist_seabed)
    factory.synchronize()
    
    
             # Ajouter un champ "Min"
    min_field_id = model.mesh.field.add("Min")

        # Définir les paramètres pour le champ "Min"
    model.mesh.field.setNumbers(min_field_id, "FieldsList", [threshold_cable,threshold_id4])
        # Définir le champ "Min" comme champ de taille de maillage
    model.mesh.field.setAsBackgroundMesh(min_field_id)

        # Synchroniser le modèle
    factory.synchronize()



    # Générer le maillage
    model.mesh.generate(2)





    # MshFileVersion settings : do not change !
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        
    # Mesh save
    meshName_msh = meshName + ".msh"
    gmsh.write(meshName_msh)
    #visua
    if visual:
        gmsh.fltk.run()
    gmsh.finalize()
    return
    



def convertMsh2Mail(meshName):    
    """

    Parameters
    ----------
    meshName : mesh and job name for the METIS computation

    Returns
    -------
    The "meshName".mail METIS mesh file 
    
    Conversion of gmsh mesh and METIS mesh

    """
    # ver
    # libs_gfortran = ['gfortran']

    import os
    script_dir = os.path.dirname(os.path.abspath(__file__))
    mshTomail_path = os.path.join(script_dir, 'mshTomail')
    mail_gmsh_2_met_path = os.path.join(script_dir, 'mail_gmsh_2_met.f90')

    if not os.path.isfile(mshTomail_path):
        subprocess.call(["gfortran", "-o", mshTomail_path, mail_gmsh_2_met_path])

    subprocess.call([
        mshTomail_path,
        meshName + ".msh",
        meshName + ".mail",
        meshName + ".geom",
        "2D/3D"
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return

