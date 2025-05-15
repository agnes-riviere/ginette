# -*- coding: utf-8 -*-
"""
This module provides utilities for performing grid search over zone parameters and modifying zone parameter files
for simulation studies.
Functions:
- generate_grid_search_table(zones_to_invert, parameters, param_ranges):
    Generates a DataFrame containing all combinations of specified parameters for each zone to be inverted,
    suitable for grid search simulations.
- modify_zone_parameters(param_structure, zones_to_invert, parameters_to_invert, value_parameter):
    Modifies the parameters of specified zones in the 'E_zone_parameter.dat' file according to provided values,
    updating only the selected parameters for each zone.
- run_one_simulation_2D(repertory, id_sim, nb_zone, parameters, value_zone_parameter):
    Runs a single 2D simulation, then renames and saves the output file with a simulation-specific name.
"""

import pandas as pd
import numpy as np
from itertools import product

def generate_grid_search_table(zones_to_invert, parameters, param_ranges):
    """
    Generates a parameter table for each simulation in a grid search.

    Parameters:
    - zones_to_invert (list): List of zone numbers to invert (e.g., [4, 5]).
    - parameters (list): List of parameter names to invert (e.g., ['k', 'n', 'l']).
    - param_ranges (dict): Dictionary containing ranges (min, max, step) for each parameter, 
                           where keys are zone numbers and values are dictionaries of parameter ranges.
                           EXAMPLE:
                           {
                               4: {'k': (1, 10, 1), 'n': (0.1, 0.5, 0.1)},
                               5: {'k': (1, 10, 1), 'n': (0.1, 0.5, 0.1)}
                           }

    Returns:
    - pd.DataFrame: Table containing all parameter combinations for each simulation.
    """
    import pandas as pd
    import numpy as np
    from itertools import product

    # Generate parameter ranges for each zone
    zone_param_values = {}
    for zone in zones_to_invert:
        zone_ranges = param_ranges.get(zone, {})
        param_values = []
        for param in parameters:
            if param in zone_ranges:
                min_val, max_val, step = zone_ranges[param]
                param_values.append(np.arange(min_val, max_val + step, step))
        if param_values:
            zone_param_values[zone] = list(product(*param_values))

    # Generate all combinations across zones
    all_combinations = list(product(*zone_param_values.values()))

    # Flatten the combinations and create a DataFrame
    flattened_combinations = []
    for combination in all_combinations:
        flattened_dict = {}
        for zone, zone_combination in zip(zones_to_invert, combination):
            for param, value in zip(parameters, zone_combination):
                flattened_dict[f"{param}{zone}"] = value
        flattened_combinations.append(flattened_dict)

    # Convert to DataFrame
    param_table = pd.DataFrame(flattened_combinations)

    return param_table

def modify_zone_parameters(param_structure, zones_to_invert, parameters_to_invert, value_parameter):
    """
    Modifie les paramètres de zones dans le fichier 'E_zone_parameter.dat' 
    en fonction des zones et paramètres spécifiés.

    Paramètres :
    - param_structure : liste des noms de colonnes (ex. ['zone','k', 'n', 'l', 'cpm', 'r'])
    - zones_to_invert : liste des numéros de zones à modifier (indices à partir de 1)
    - parameters_to_invert : liste des paramètres à modifier (ex. ['k', 'n'])
    - value_parameter : Series contenant les nouvelles valeurs, avec des colonnes comme 'k4', 'n4'
    """

    # Lire les lignes du fichier
    with open("E_zone_parameter.dat", "r") as f:
        lines = f.read().splitlines()

    # Parcourir les zones à modifier
    for zone in zones_to_invert:
        line_index = zone - 1  # index Python (0-based)

        if line_index >= len(lines):
            print(f" Zone {zone} absente du fichier.")
            continue

        columns = lines[line_index].split()

        for param in parameters_to_invert:
            if param not in param_structure:
                print(f" Paramètre '{param}' absent de la structure.")
                continue

            param_index = param_structure.index(param)
            col_name = f"{param}{zone}"

            if col_name not in value_parameter.index:
                print(f" Colonne '{col_name}' absente dans la Series de valeurs.")
                continue

            value = value_parameter[col_name]  # Pas de .values[0] ici
            formatted_value = f"{value:.8e}".replace('e', 'd')  # Format FORTRAN
            columns[param_index] = formatted_value

        lines[line_index] = '\t'.join(columns)

    # Réécriture du fichier avec les lignes modifiées
    with open("E_zone_parameter.dat", "w") as f:
        f.write('\n'.join(lines) + '\n')

def run_one_simulation_2D(repertory,id_sim,nb_zone,parameters,value_zone_parameter):
    run_direct_model_2D()
    # save the file S_temp_PT100_t.dat with the name of the simulation at the beginning and _
    # Construct the new file name
    new_file_name = f"{repertory}/{id_sim}_S_temp_PT100_t.dat"
    
    # Save the file with the new name
    os.rename("S_temp_PT100_t.dat", new_file_name)
