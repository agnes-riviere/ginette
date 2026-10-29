cmake_minimum_required(VERSION 3.10)
project(GinetteProject Fortran)

# Ajoutez vos fichiers source Fortran ici
set(SOURCES
    ../../src/ginette_V2.f90
)

# Créez l'exécutable
add_executable(ginette ${SOURCES})
