cmake_minimum_required(VERSION 3.12)
project(ginette_project LANGUAGES Fortran)

# Active les options d’optimisation
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")

# Chemin vers ton code source
set(SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../src/ginette_V2.f90")

# Définis l’exécutable à générer
add_executable(ginette ${SRC_PATH})

# Dossier de sortie de l’exécutable 
set_target_properties(ginette PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/input_ginette"
)
