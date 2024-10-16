#!/usr/bin/env python


import os
import glob

def prepare_ginette_directories(base_path, subdirectories=['SENSI', 'OUTPUT']):
    """
    Change le répertoire de travail pour le modèle Ginette et crée les répertoires de sortie nécessaires.

    Args:
        base_path (str): Le chemin vers le répertoire contenant le modèle Ginette.
        subdirectories (list): Une liste de sous-répertoires à créer dans le répertoire de base.
    """
    os.chdir(base_path)
    print("Current working directory: {0}".format(os.getcwd()))

    for subdir in subdirectories:
        if not os.path.exists(subdir):
            os.makedirs(subdir)
            print(f"Directory '{subdir}' created.")
        else:
            print(f"Directory '{subdir}' already exists.")
            
    # if files exist in OUTPUT or SENSI delete them
    files = glob.glob('OUTPUT/*')
    for f in files:
        os.remove(f)
        print("file deleted",f)
    files = glob.glob('SENSI/*')
    for f in files:
        os.remove(f)
        print("file deleted",f)        
            
def compile_ginette():
    """
    Compile Ginette si le fichier exécutable n'existe pas.
    """
    if os.path.isfile('ginette'):
        print("ginette exists")
    else:
        print("ginette does not exist")
        print("You must compile ginette in the current directory")
        print("Command: gfortran -o ginette ../../src/ginette_V2.f")
        os.system("gfortran -o ginette ../../src/ginette_V2.f")  # Utilisation de os.system
