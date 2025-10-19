#!/bin/bash
# filepath: /home/ariviere/Programmes/ginette/application/2024_TD_ENS/GINETTE_SENSI/compile.sh

echo "=== Compilation de Ginette ==="

# Méthode 1: Compilation avec CMake
if [ -f "CMakeLists.txt" ]; then
    echo "Compilation avec CMake..."
    
    # Créer le répertoire build s'il n'existe pas
    mkdir -p build
    cd build
    
    # Configuration et compilation
    cmake .. && cmake --build .
    
    # Retourner au répertoire principal
    cd ..
    
    if [ -f "ginette" ]; then
        echo "✓ Compilation réussie avec CMake"
        exit 0
    fi
fi

# Méthode 2: Compilation directe avec gfortran
echo "Compilation directe avec gfortran..."

# Trouver tous les fichiers Fortran
FORTRAN_FILES=$(find . -name "*.f90" -o -name "*.f" -o -name "*.F90" -o -name "*.F" | grep -v build)

if [ -n "$FORTRAN_FILES" ]; then
    echo "Fichiers Fortran trouvés:"
    echo "$FORTRAN_FILES"
    
    # Compilation
    gfortran -O2 -fdefault-real-8 -fdefault-double-8 -o ginette $FORTRAN_FILES
    
    if [ -f "ginette" ]; then
        echo "✓ Compilation réussie avec gfortran"
        exit 0
    fi
else
    echo "❌ Aucun fichier Fortran trouvé"
fi

echo "❌ Échec de la compilation"
exit 1