# Creation date: 2020-07-15
#auteur: Agnes Riviere
import subprocess
import os
import glob
import platform
from pathlib import Path


def ginette_executable_name():
    """Nom de l'exécutable Ginette selon la plateforme (ginette.exe sous Windows)."""
    return 'ginette.exe' if platform.system() == 'Windows' else 'ginette'


def setup_lab_environment():
    """
    Prépare l'environnement du TD : localise le dossier SYNTHETIC_CASES à côté du
    notebook, y prépare le dossier de travail, et compile l'exécutable 'ginette'
    s'il n'existe pas déjà.

    À appeler juste après avoir ajouté 'src_python' au sys.path.
    Peut être appelée plusieurs fois de suite dans le même kernel (idempotent) :
    comme elle change le répertoire de travail courant vers SYNTHETIC_CASES,
    un ré-appel ne duplique pas ce chemin.

    Returns:
        (Path, Path): (dossier du TD, dossier SYNTHETIC_CASES)
    """
    cwd = Path.cwd().resolve()
    td_dir = cwd.parent if cwd.name == 'SYNTHETIC_CASES' else cwd
    synthetic_cases_dir = td_dir / 'SYNTHETIC_CASES'

    prepare_ginette_directories(str(synthetic_cases_dir))
    compile_ginette()

    return td_dir, synthetic_cases_dir

def prepare_ginette_directories(base_path, subdirectories=None):
    """
    Change le répertoire de travail pour le modèle Ginette et crée les éventuels
    répertoires de sortie nécessaires.

    Args:
        base_path (str): Le chemin vers le répertoire contenant le modèle Ginette.
        subdirectories (list): Sous-répertoires à créer/nettoyer dans le répertoire de base.
            Aucun par défaut : Ginette écrit ses fichiers de sortie (Sim_*.dat) directement
            dans base_path, SENSI/ et OUTPUT/ ne sont pas utilisés par ce workflow.
    """
    os.chdir(base_path)
    print("Répertoire de travail courant : {0}".format(os.getcwd()))

    for subdir in (subdirectories or []):
        if not os.path.exists(subdir):
            os.makedirs(subdir)
            print(f"Répertoire '{subdir}' créé.")
        else:
            print(f"Répertoire '{subdir}' déjà existant.")

        files = glob.glob(f'{subdir}/*')
        for f in files:
            os.remove(f)
            print("fichier supprimé", f)


def _find_fortran_sources():
    """
    Cherche les fichiers source Fortran dans les répertoires candidats habituels.
    Retourne (fortran_files, source_dir), ([], None) si rien n'est trouvé.
    """
    possible_source_dirs = [
        '.',              # Répertoire courant
        'src',            # Sous-répertoire src
        '../src',         # Répertoire src parent
        '../../src',      # Répertoire src grand-parent
        '../../../src',   # Répertoire src (racine du dépôt depuis SYNTHETIC_CASES)
        'source',         # Sous-répertoire source
        'fortran',        # Sous-répertoire fortran
    ]

    for search_dir in possible_source_dirs:
        if os.path.exists(search_dir):
            temp_files = []
            for ext in ['*.f90', '*.f', '*.F90', '*.F']:
                temp_files.extend(glob.glob(os.path.join(search_dir, ext)))

            if temp_files:
                return temp_files, search_dir
    return [], None


def _exe_is_stale(exe_name, fortran_files):
    """
    True si l'exécutable n'existe pas, ou si un des fichiers source est plus
    récent que lui. Sans ce test, un exécutable déjà présent (session
    précédente, `git pull` qui met à jour le Fortran sans toucher au binaire)
    n'est jamais recompilé : le TD tourne alors avec une version figée du
    code, silencieusement.
    """
    if not os.path.isfile(exe_name):
        return True
    exe_mtime = os.path.getmtime(exe_name)
    return any(os.path.getmtime(f) > exe_mtime for f in fortran_files)


def compile_ginette():
    """
    Compile Ginette si le fichier exécutable n'existe pas ou est plus ancien
    que le source Fortran (voir _exe_is_stale).
    Fonctionne sur Linux, macOS et Windows (l'extension .exe est gérée automatiquement).
    """
    exe_name = ginette_executable_name()
    fortran_files, source_dir = _find_fortran_sources()

    if os.path.isfile(exe_name) and not _exe_is_stale(exe_name, fortran_files):
        print(f"{exe_name} existe déjà et est à jour")
    else:
        print(f"{exe_name} n'existe pas ou est obsolète")
        print("Tentative de compilation de ginette...")

        if not fortran_files:
            if os.path.isfile(exe_name):
                print("Aucun fichier source Fortran trouvé : conservation de l'exécutable existant.")
                return
            print("Aucun fichier source Fortran trouvé dans les répertoires attendus")
            print("Vérifiez que les fichiers source Fortran (.f90, .f, .F90, .F) sont présents")
            return
        print(f"Fichiers Fortran trouvés dans : {source_dir}")

        # Méthode 1: Essayer avec le script de compilation (Linux/macOS uniquement)
        if platform.system() != 'Windows' and os.path.isfile('compile.sh'):
            try:
                print("Found compile.sh, attempting compilation...")
                subprocess.run(['chmod', '+x', 'compile.sh'], check=True)
                result = subprocess.run(['./compile.sh'], capture_output=True, text=True, timeout=300)
                if result.returncode == 0 and os.path.isfile(exe_name):
                    print("✓ Ginette compiled successfully with compile.sh")
                    return
                else:
                    print(f"compile.sh failed: {result.stderr}")
            except Exception as e:
                print(f"compile.sh execution failed: {e}")
        
        # Méthode 2: CMake avec création automatique du CMakeLists.txt si nécessaire
        if not os.path.isfile('CMakeLists.txt') and fortran_files:
            print("Creating CMakeLists.txt automatically...")
            cmake_content = f'''cmake_minimum_required(VERSION 3.10)
project(ginette LANGUAGES Fortran)

set(CMAKE_Fortran_FLAGS "${{CMAKE_Fortran_FLAGS}} -O2 -fdefault-real-8 -fdefault-double-8")

set(FORTRAN_SOURCES
'''
            for f in fortran_files:
                cmake_content += f'    "{f}"\n'
            
            cmake_content += ''')

add_executable(ginette ${FORTRAN_SOURCES})
set_target_properties(ginette PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
'''
            
            with open('CMakeLists.txt', 'w') as f:
                f.write(cmake_content)
        
        if os.path.isfile('CMakeLists.txt'):
            try:
                print("Found CMakeLists.txt, attempting CMake compilation...")
                if not os.path.exists('build'):
                    os.makedirs('build')
                
                result1 = subprocess.run(['cmake', '-B', 'build', '.'],
                                       capture_output=True, text=True, timeout=120)
                if result1.returncode == 0:
                    result2 = subprocess.run(['cmake', '--build', 'build'],
                                           capture_output=True, text=True, timeout=300)
                    if result2.returncode == 0 and os.path.isfile(exe_name):
                        print("✓ Ginette compiled successfully with CMake")
                        return
                    else:
                        print(f"CMake build failed: {result2.stderr}")
                else:
                    print(f"CMake configuration failed: {result1.stderr}")
            except Exception as e:
                print(f"CMake compilation failed: {e}")
        
        # Méthode 3: gfortran direct
        try:
            print(f"Attempting direct gfortran compilation with files from {source_dir}...")
            print(f"Compiling files: {[os.path.basename(f) for f in fortran_files]}")

            cmd = ['gfortran', '-O2', '-fdefault-real-8', '-fdefault-double-8',
                   '-o', exe_name] + fortran_files
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if result.returncode == 0 and os.path.isfile(exe_name):
                print("✓ Ginette compiled successfully with gfortran")
                return
            else:
                print(f"Gfortran compilation failed: {result.stderr}")
        except Exception as e:
            print(f"Gfortran compilation failed: {e}")

        # Instructions manuelles
        print("\n❌ Automatic compilation failed.")
        print("Manual compilation options:")
        if fortran_files:
            print(f"1. Direct compilation: gfortran -O2 -o {exe_name} {' '.join(fortran_files)}")
        print("2. Check if there's a specific Makefile or build instructions")
        print(f"Current directory: {os.getcwd()}")
        print(f"Found source files in '{source_dir}': {[os.path.basename(f) for f in fortran_files]}")