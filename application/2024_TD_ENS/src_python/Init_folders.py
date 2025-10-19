# Creation date: 2020-07-15
#auteur: Agnes Riviere
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
    """
    Compile Ginette si le fichier exécutable n'existe pas.
    """
    if os.path.isfile('ginette'):
        print("ginette exists")
    else:
        print("ginette does not exist")
        print("Attempting to compile ginette...")
        
        # Chercher les fichiers Fortran dans différents répertoires
        possible_source_dirs = [
            '.',           # Répertoire courant
            'src',         # Sous-répertoire src
            '../src',      # Répertoire src parent
            'source',      # Sous-répertoire source
            'fortran',     # Sous-répertoire fortran
        ]
        
        fortran_files = []
        source_dir = None
        
        for search_dir in possible_source_dirs:
            if os.path.exists(search_dir):
                temp_files = []
                for ext in ['*.f90', '*.f', '*.F90', '*.F']:
                    temp_files.extend(glob.glob(os.path.join(search_dir, ext)))
                
                if temp_files:
                    fortran_files = temp_files
                    source_dir = search_dir
                    print(f"Found Fortran files in: {search_dir}")
                    break
        
        if not fortran_files:
            print("No Fortran source files found in any expected directory")
            print(f"Searched in: {possible_source_dirs}")
            print("Please ensure Fortran source files (.f90, .f, .F90, .F) are present")
            return
        
        # Méthode 1: Essayer avec le script de compilation
        if os.path.isfile('compile.sh'):
            try:
                print("Found compile.sh, attempting compilation...")
                subprocess.run(['chmod', '+x', 'compile.sh'], check=True)
                result = subprocess.run(['./compile.sh'], capture_output=True, text=True, timeout=300)
                if result.returncode == 0 and os.path.isfile('ginette'):
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
                    if result2.returncode == 0 and os.path.isfile('ginette'):
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
                   '-o', 'ginette'] + fortran_files
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            if result.returncode == 0:
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
            print(f"1. Direct compilation: gfortran -O2 -o ginette {' '.join(fortran_files)}")
        print("2. Check if there's a specific Makefile or build instructions")
        print(f"Current directory: {os.getcwd()}")
        print(f"Found source files in '{source_dir}': {[os.path.basename(f) for f in fortran_files]}")