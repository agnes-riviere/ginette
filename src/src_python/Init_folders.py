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
_PE_MAGIC = b'MZ'


def _ginette_binary_matches_host(path='ginette'):
    """
    True si le binaire ginette existe et que ses premiers octets
    correspondent au format natif de la plateforme hôte - ne l'exécute pas,
    se contente de lire son en-tête :
    - macOS (Darwin) : Mach-O
    - Windows : PE ("MZ") - la compilation gfortran produit en réalité
      "ginette.exe" (pas "ginette") sur Windows, donc ce chemin ne sert que
      si un ginette.exe a été renommé sans l'extension ; non testé en
      pratique, ce projet n'étant pas utilisé sous Windows à notre
      connaissance.
    - tout le reste (Linux et autres Unix - BSD, etc.) : ELF, format
      produit par gfortran sur toutes ces plateformes.
    """
    if not os.path.isfile(path):
        return False
    try:
        with open(path, 'rb') as f:
            magic = f.read(4)
    except OSError:
        return False
    system = platform.system()
    if system == 'Darwin':
        return magic in _MACHO_MAGICS
    if system == 'Windows':
        return magic[:2] == _PE_MAGIC
    return magic == _ELF_MAGIC  # Linux et autres Unix


def _ginette_binary_is_fresh(path, source_path):
    """
    True si le binaire ginette correspond à la plateforme hôte ET a été
    compilé APRES la dernière modification du source Fortran. Sans ce
    deuxième test, un exécutable laissé sur un poste ou un cluster (session
    précédente, `git pull` qui met à jour ginette_V2.f90 sans toucher au
    binaire déjà présent) n'est jamais recompilé : on tourne alors avec une
    version figée du code sans le savoir, ce qui peut ressembler à du bruit
    numérique alors que c'est une simple exécution obsolète.
    """
    if not _ginette_binary_matches_host(path):
        return False
    try:
        return os.path.getmtime(path) >= os.path.getmtime(source_path)
    except OSError:
        return False  # source introuvable : on laisse (re)compiler échouer plus loin avec un message clair


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
    (Re)compile Ginette si l'exécutable est absent, incompatible avec la
    plateforme hôte, ou plus ancien que ginette_V2.f90 (voir
    _ginette_binary_is_fresh).
    """
    source = '../../src/ginette_V2.f90'
    if _ginette_binary_is_fresh('ginette', source):
        print("ginette exists")
    else:
        print("ginette does not exist, does not match this platform, or is older than the source - (re)compiling")
        result = subprocess.run(['gfortran', '-o', 'ginette', source])
        if result.returncode != 0:
            raise RuntimeError(f"gfortran compilation of ginette failed (exit code {result.returncode})")
        print("ginette compiled")


def compile_ginette_src(dir_ginette, flags=()):
    """
    (Re)compile Ginette si l'exécutable est absent, incompatible avec la
    plateforme hôte - typiquement un binaire Linux ELF copié tel quel sur
    macOS, qui échoue silencieusement en "Exec format error" à l'exécution -
    ou plus ancien que ginette_V2.f90 (voir _ginette_binary_is_fresh).

    flags : options supplémentaires passées à gfortran (ex. ("-O2",) pour les
    cas 2D longs). Sans effet si l'exécutable existe déjà et est à jour : le
    supprimer pour forcer la recompilation avec d'autres flags.
    """
    source = dir_ginette + '/src/ginette_V2.f90'
    if _ginette_binary_is_fresh('ginette', source):
        print("ginette exists")
    else:
        print("ginette does not exist, does not match this platform, or is older than the source - (re)compiling")
        result = subprocess.run(['gfortran', *flags, '-o', 'ginette', source])
        if result.returncode != 0:
            raise RuntimeError(f"gfortran compilation of ginette failed (exit code {result.returncode})")
        print("ginette compiled")

