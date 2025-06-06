# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 14:57:23 2023
Metro.py
This module provides utility functions for handling metro simulation parameters, 
reading input files, computing marginal errors, and calculating vector norms. 
It is designed to work with the Ginette direct model for metro simulations.
Functions:
----------
- calcul_marginale(df, parametre, erreur_colonne='Erreur'):
    Groups a DataFrame by a specified parameter and sums the error column to compute marginal errors.
- norme_2(vect1, vect2, sigma):
    Computes the normalized sum of squared differences (L2 norm) between two vectors.
- read_input_metro():
    Reads simulation configuration from 'read_input_metro.txt', including number of iterations, 
    measurement uncertainties, number of zones, and parameter names.
- read_bounds_params(nb_zone, parameters):
    Reads parameter bounds and standard deviations for each zone from 'bound_params.txt', 
    returning a list of dictionaries mapping parameters to their bounds and std.
Dependencies:
-------------
- numpy
- pandas
- os
- src.src_python.Direct_model (setup_ginette, generate_zone_parameters, initial_conditions, boundary_conditions, run_direct_model, remove_first_two_days)
"""
import numpy as np
import pandas as pd
import os
from src.src_python.Direct_model import setup_ginette, generate_zone_parameters, initial_conditions, boundary_conditions, run_direct_model, remove_first_two_days



def calcul_marginale(df, parametre, erreur_colonne='Erreur'):
    # Regrouper par le paramètre choisi et sommer les erreurs
    marginale = (
        df.groupby(parametre)[erreur_colonne]
        .sum()
        .reset_index()
        .sort_values(by=parametre)  # Trier par ordre croissant
    )
    return marginale




def norme_2(vect1,vect2,sigma):
    """
    Calcule la somme des carrés des différences entre deux vecteurs.
    
    Parameters:
        vect1 (list ou numpy.array) : Premier vecteur
        vect2 (list ou numpy.array) : Deuxième vecteur
        
    Returns:
        float : La somme des carrés des différences
    """
    if len(vect1) != len(vect2):
        raise ValueError("Les deux vecteurs doivent avoir la même taille.")
    
    return (sum((a - b) ** 2 for a, b in zip(vect1, vect2)))/sigma**2



def read_input_metro():
    """
    Reads the file 'read_input_metro.txt'.

    The file should contain:
    - First line: Number of iterations (integer)
    - Second line: Uncertainties of the temperature measurements (space-separated floats)
    - Third line: Number of zones (integer)
    - Fourth line: List of parameter names (space-separated strings, e.g., ['k', 'n', 'lambda', 'c'])

    Returns:
    - iterations (int): Number of iterations
    - uncertainties (list of float): List of uncertainties
    - nb_zone (int): Number of zones
    - parameters (list of str): List of parameter names
    """
    try:
        with open("read_input_metro.txt", 'r') as f:
            lines = f.readlines()
            if len(lines) < 4:
                raise ValueError("Error: 'read_input_metro.txt' must contain at least 4 lines.")
            iterations = int(lines[0].strip())
            uncertainties = list(map(float, lines[1].strip().split()))
            nb_zone = int(lines[2].strip())
            parameters = lines[3].strip().split()
    except FileNotFoundError:
        raise FileNotFoundError("Error: 'read_input_metro.txt' file not found.")
    except ValueError as e:
        raise ValueError(f"Error processing file: {e}")

    return iterations, uncertainties, nb_zone, parameters




def read_bounds_params(nb_zone, parameters):
    """
    Reads the file 'bound_params.txt' containing zone, parameter, min_value, max_value, and std.

    Parameters:
    - nb_zone: Number of zones
    - parameters: List of parameter names (e.g., ['k', 'n', 'lambda', 'c'])

    Returns:
    - A list of dictionaries, one for each zone, where each dictionary maps parameters to their bounds and std.
    """
    params_dict = {zone: {param: None for param in parameters} for zone in range(1, nb_zone + 1)}

    try:
        with open("bound_params.txt", 'r') as f:
            for line in f:
                if not line.strip():
                    continue

                columns = line.strip().split()
                if len(columns) != 5:
                    print(f"Invalid line format: {line}")
                    continue

                zone, param, min_value, max_value, std = int(columns[0]), columns[1], float(columns[2]), float(columns[3]), float(columns[4])

                if 1 <= zone <= nb_zone and param in parameters:
                    params_dict[zone][param] = (min_value, max_value, std)
                else:
                    print(f"Invalid entry: {line}")

    except FileNotFoundError:
        print("Error: 'bound_params.txt' file not found.")
    except ValueError as e:
        print(f"Error processing line: {line}. {e}")

    return [params_dict[zone] for zone in range(1, nb_zone + 1)]


