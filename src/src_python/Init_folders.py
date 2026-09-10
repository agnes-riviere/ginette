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
import platform

# Magic numbers (premiers octets du fichier) identifiant le format natif de
# chaque plateforme - sert à détecter un binaire ginette compilé pour une
# AUTRE plateforme (typiquement Linux ELF copié tel quel sur macOS) avant de
# tenter de l'exécuter, ce qui échoue silencieusement en "Exec format error"
# plutôt que de déclencher la recompilation attendue.
_ELF_MAGIC = b'\x7fELF'
_MACHO_MAGICS = (b'\xcf\xfa\xed\xfe', b'\xce\xfa\xed\xfe', b'\xca\xfe\xba\xbe', b'\xbe\xba\xfe\xca')


def _ginette_binary_matches_host(path='ginette'):
    """
    True si le binaire ginette existe et que ses premiers octets
    correspondent au format natif de la plateforme hôte (ELF sur Linux,
    Mach-O sur macOS) - ne l'exécute pas, se contente de lire son en-tête.
    """
    if not os.path.isfile(path):
        return False
    try:
        with open(path, 'rb') as f:
            magic = f.read(4)
    except OSError:
        return False
    system = platform.system()
    if system == 'Linux':
        return magic == _ELF_MAGIC
    if system == 'Darwin':
        return magic in _MACHO_MAGICS
    return True  # plateforme non reconnue : on ne bloque pas, laisse échouer à l'exécution

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
    (Re)compile Ginette si l'exécutable est absent OU incompatible avec la
    plateforme hôte (voir _ginette_binary_matches_host).
    """
    if _ginette_binary_matches_host('ginette'):
        print("ginette exists")
    else:
        print("ginette does not exist or does not match this platform - (re)compiling")
        subprocess.run(['gfortran', '-o', 'ginette', '../../src/ginette_V2.f90'])
        if os.path.isfile('ginette'):
            print("ginette compiled")


def compile_ginette_src(dir_ginette):
    """
    (Re)compile Ginette si l'exécutable est absent OU incompatible avec la
    plateforme hôte - typiquement un binaire Linux ELF copié tel quel sur
    macOS, qui échoue silencieusement en "Exec format error" à l'exécution
    si on se contente de tester sa présence (voir _ginette_binary_matches_host).
    """
    if _ginette_binary_matches_host('ginette'):
        print("ginette exists")
    else:
        print("ginette does not exist or does not match this platform - (re)compiling")
        subprocess.run(['gfortran', '-o', 'ginette', dir_ginette + '/src/ginette_V2.f90'])
        if os.path.isfile('ginette'):
            print("ginette compiled")

