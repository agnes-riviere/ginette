cmake_minimum_required(VERSION 3.10)
project(ginette LANGUAGES Fortran)

# Set the Fortran compiler flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2 -fdefault-real-8 -fdefault-double-8")

# Specify the exact source file location
set(GINETTE_SOURCE "/home/ariviere/Programmes/ginette/src/ginette_V2.f90")

# Check if the source file exists
if(NOT EXISTS ${GINETTE_SOURCE})
    message(FATAL_ERROR "Source file not found: ${GINETTE_SOURCE}")
endif()

# Also check for any additional Fortran files in the same directory
get_filename_component(SOURCE_DIR ${GINETTE_SOURCE} DIRECTORY)
file(GLOB ADDITIONAL_SOURCES 
    "${SOURCE_DIR}/*.f90" 
    "${SOURCE_DIR}/*.f" 
    "${SOURCE_DIR}/*.F90" 
    "${SOURCE_DIR}/*.F"
)

# Remove duplicates and ensure main source is included
list(REMOVE_DUPLICATES ADDITIONAL_SOURCES)
set(FORTRAN_SOURCES ${GINETTE_SOURCE})

# Add any additional sources if found (excluding the main one to avoid duplicates)
foreach(source ${ADDITIONAL_SOURCES})
    if(NOT ${source} STREQUAL ${GINETTE_SOURCE})
        list(APPEND FORTRAN_SOURCES ${source})
    endif()
endforeach()

# Create the executable
add_executable(ginette ${FORTRAN_SOURCES})

# Set output directory to current directory
set_target_properties(ginette PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# Print information about found files
message(STATUS "Main Ginette source: ${GINETTE_SOURCE}")
message(STATUS "All Fortran source files:")
foreach(file ${FORTRAN_SOURCES})
    message(STATUS "  ${file}")
endforeach()
