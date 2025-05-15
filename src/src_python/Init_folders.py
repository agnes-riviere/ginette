# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 14:23:45 2023
This module provides utility functions to prepare the working directories and compile the Ginette model.
Functions:
    prepare_ginette_directories(base_path, subdirectories=['SENSI', 'OUTPUT']):
        Changes the working directory to the specified base path, creates necessary subdirectories
        ('SENSI' and 'OUTPUT' by default), and cleans up any existing files in these directories.
    compile_ginette():
        Compiles the Ginette Fortran executable using gfortran if it does not already exist in the current directory.
"""
import subprocess
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
    import shutil
    import subprocess
    """
    Compile Ginette if the executable does not exist.
    This function uses Meson for building the project.
    """
    if os.path.isfile('ginette'):
        print("ginette exists")
    else:
        print("ginette does not exist")
        subprocess.run(['gfortran', '-o', 'ginette', '../../src/ginette_V2.f90'])
        if os.path.isfile('ginette'):
            print("ginette compiled")