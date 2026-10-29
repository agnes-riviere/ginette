# -*- coding: utf-8 -*-
"""
This module provides utilities for performing grid search over zone parameters and modifying zone parameter files
for simulation studies.
Main Functionalities:
---------------------
1. Grid Search Table Generation:
    - `generate_grid_search_table(zones_to_invert, parameters, param_ranges)`:
      Generates a pandas DataFrame containing all possible combinations of specified parameters for each zone,
      suitable for exhaustive grid search in simulation studies.
2. Zone Parameter Modification:
    - `modify_zone_parameters(param_structure, zones_to_invert, parameters_to_invert, value_parameter)`:
      updating only the selected parameters for each zone. Handles log-to-linear conversion for 'k' parameters.
3. Simulation Execution:
    - `run_one_simulation_2D(...)`:
      Runs a single 2D simulation for a given station and parameter set, computes the mean squared error (MSE)
      between simulated and observed temperatures, and saves simulation results for further analysis.
4. Marginal and Joint Posterior Visualization:
    - `plot_marginal_1D(df, param_name, param_units=None)`:
      Plots the 1D marginal posterior distribution for a given parameter, with optional unit annotation.
    - `plot_marginal_2D(df, param_x, param_y, param_units=None, cmap="viridis")`:
      Plots the 2D marginal posterior distribution for a pair of parameters as a heatmap, with optional units.
    - `plot_joint_posterior(df, param_cols, param_units=None, cmap="viridis")`:
      Creates a grid of plots showing 1D and 2D marginal posteriors for a set of parameters, enabling visualization
      of parameter interactions and most probable values.
Dependencies:
-------------
- numpy
- pandas
- matplotlib
- seaborn
- fortranformat
- Custom modules: Init_folders, Direct_model, Read_obs, Plot, stat_critere, mesh_generator
Usage Example:
--------------
Typical workflow involves:
1. Generating a grid search table for parameter combinations.
2. Iteratively modifying zone parameters and running simulations.
3. Collecting simulation results and computing statistical criteria (e.g., MSE).
4. Visualizing marginal and joint posterior distributions to analyze parameter sensitivities and uncertainties.
------
- The module assumes the existence of specific data files and directory structures.
- Functions are designed for extensibility and integration with larger simulation pipelines.
- Visualization functions are intended for exploratory data analysis and result interpretation.

"""
import os
import importlib
from itertools import product

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import fortranformat as ft

from src.src_python import Init_folders
from src.src_python import Direct_model
from src.src_python import Read_obs
from src.src_python import Plot
from src.src_python import stat_critere


# Import all functions/classes from the relevant modules
from src.src_gmsh.mesh_generator import *
from src.src_python.Init_folders import *
from src.src_python.Direct_model import *
from src.src_python.Read_obs import *
from src.src_python.Plot import *
from src.src_python.stat_critere import *




importlib.reload(Init_folders)
importlib.reload(Direct_model)


def generate_grid_search_table(zones_to_invert, parameters, param_ranges):
    """
    Generates a DataFrame containing all possible combinations of specified parameters for each zone,
    suitable for grid search in simulation studies.

    Parameters:
        zones_to_invert (list of int): 
            List of zone numbers to invert (e.g., [4, 5]). Each zone will have its own set of parameter values.
        parameters (list of str): 
            List of parameter names to vary (e.g., ['k', 'n', 'l']). These parameters will be combined for each zone.
        param_ranges (dict): 
            Dictionary specifying the range or explicit values for each parameter in each zone.
            The keys are zone numbers (int), and the values are dictionaries mapping parameter names (str)
            to either:
                - a tuple (min, max, step) for generating values using np.arange, or
                - a list/array of explicit values.

    Returns:
        pd.DataFrame: 
            DataFrame where each row represents a unique simulation configuration.
            Columns are named as 'parameterNameZoneNumber' (e.g., 'k4', 'n4', 'k5', 'n5'),
            containing all possible combinations of parameter values across the specified zones.
    """

    # Prepare a dictionary to hold all possible value combinations for each zone
    zone_param_values = {}
    for zone in zones_to_invert:
        zone_ranges = param_ranges.get(zone, {})
        param_values = []
        # For each parameter, generate the list of values it can take for this zone
        for param in parameters:
            if param in zone_ranges:
                range_val = zone_ranges[param]
                if isinstance(range_val, (tuple, list)) and len(range_val) == 3:
                    min_val, max_val, step = range_val
                    # Compute the number of steps so that spacing is always constant and both ends are included
                    num_steps = int(round((max_val - min_val) / step)) + 1
                    values = np.linspace(min_val, max_val, num_steps)
                else:
                    values = np.array(range_val)
                # Round to 4 decimals and convert -0.0 to 0.0
                values = np.round(values, 4)
                values = np.where(np.isclose(values, 0), 0.0, values)
                param_values.append(values)
                print(f"Zone {zone}, Parameter '{param}': Values = {values}")
        if param_values:
            # Compute the cartesian product of all parameter values for this zone
            zone_param_values[zone] = list(product(*param_values))

    # Compute the cartesian product across all zones to get every possible simulation configuration
    all_combinations = list(product(*zone_param_values.values()))

    # Flatten the nested tuples into dictionaries with keys like 'k4', 'n4', etc.
    flattened_combinations = []
    for combination in all_combinations:
        flattened_dict = {}
        for zone, zone_combination in zip(zones_to_invert, combination):
            for param, value in zip(parameters, zone_combination):
                # Create column names by combining parameter and zone number
                flattened_dict[f"{param}{zone}"] = value
        flattened_combinations.append(flattened_dict)

    # Convert the list of dictionaries to a pandas DataFrame
    param_table = pd.DataFrame(flattened_combinations)
    # Round all float values to 2 decimals for output
    param_table = param_table.round(2)
    # Save to file
    param_table.to_csv("param_table.txt", sep=",", index=False, float_format="%.2f")
    return param_table

def modify_zone_parameters(param_structure, zones_to_invert, parameters_to_invert, value_parameter):
    """
    Modifies zone parameters in the 'E_zone_parameter.dat' file based on specified zones and parameters.

    This function reads the 'E_zone_parameter.dat' file, updates the specified parameters for the given zones
    using the provided values, and writes the modified data back to the file. For parameters starting with 'k',
    the function converts their values from log10 to linear scale before updating.

    Parameters:
        -----------
        param_structure (list of str): List of parameter names defining the structure/order of parameters in each row.
        zones_to_invert (list of int): List of zone numbers for which parameters should be modified.
        parameters_to_invert (list of str): List of parameter names to be updated for the specified zones.
        value_parameter (pandas.Series): Series containing new parameter values, indexed by parameter name and zone
            (e.g., 'k1', 'alpha2', etc.). For 'k' parameters, values should be in log10 and will be converted to linear.

    Notes:
        - The function expects the file 'E_zone_parameter.dat' to exist in the current working directory.
        - The file is expected to be formatted according to the Fortran format string constructed from param_structure.
        - If a specified parameter value is missing in value_parameter, a warning is printed and the value is not updated.
        - Lines that cannot be parsed according to the format are skipped with an error message.

    Raises:
        None directly, but prints warnings and errors encountered during file reading and parameter updating.

    """
    filename="E_zone_parameter_bck.dat"
    # Log10 to linear for 'k' parameters
    for zone in zones_to_invert:
        col_name = f'k{zone}'
        if col_name in value_parameter.index:
            value = value_parameter[col_name]
            if isinstance(value, (int, float, np.number)) and np.isfinite(value):
                value_parameter.loc[col_name] = 10 ** value

    # Build format string
    nb_values = len(param_structure)-1
    fmt = "(I2," + ",".join(["1X,1PE11.2"] * nb_values) + ")"
    reader = ft.FortranRecordReader(fmt)
    writer = ft.FortranRecordWriter(fmt)

    # Read lines
    lines = []
    with open(filename, "r") as f:
        for line in f:
            if not line.strip():
                continue
            try:
                parsed = reader.read(line)
                lines.append(parsed)
            except Exception as e:
                print(f"Erreur de lecture à la ligne: {line.strip()} -> {e}")
                continue

    # Modify parameters
    for row in lines:
        zone = int(row[0])
        if zone in zones_to_invert:
            for param in parameters_to_invert:
                if param in param_structure:
                    col_idx = param_structure.index(param)
                    key = f"{param}{zone}"
                    if key in value_parameter.index:
                        row[col_idx] = value_parameter[key]
                    else:
                        print(f"Warning: '{key}' not in value_parameter")



    # Write modified data
    with open("E_zone_parameter.dat", "w") as f:
        for row in lines:
            f.write(writer.write(row).rstrip() + "\n")
        f.write(" \n")


def run_one_simulation_2D(Station,repertory,date_simul_bg,name_sensor,param_table,
                          id_sim,param_struct,zones_to_invert,parameters_to_invert,
                          param_table_simul,dir_ginette):
    """
    Runs a single 2D simulation for a given station and parameter set, computes the mean squared error (MSE)
    between simulated and observed temperatures, and saves simulation results for further analysis.

    Parameters:
        Station (str): Identifier for the station being simulated.
        repertory (str): Path to the directory containing simulation files and outputs.
        date_simul_bg (str): Simulation start date in 'YYYY-MM-DD' format.
        obs_temp (pd.DataFrame): DataFrame with observed temperature data. Columns should match sensor names.
        name_sensor (list of str): List of sensor names to use for simulation and comparison.
        param_table (pd.DataFrame): DataFrame containing all parameter sets, including an 'index_sim' column.
        id_sim (int): Simulation ID, used to select the parameter set from param_table.
        param_struct (list): List of parameter names defining the model structure.
        zones_to_invert (list): List of zone identifiers for which parameters will be inverted.
        parameters_to_invert (list): List of parameter names to be inverted during the simulation.
        param_table_simul (pd.DataFrame): DataFrame to which the current simulation's parameters will be appended.
        sigma (float): Standard deviation used in the MSE calculation.

    Returns:
        param_table_simul (pd.DataFrame): Updated DataFrame including the parameters used for this simulation.

    Notes:
        - The function modifies zone parameters using the provided structure and inversion lists.
        - It runs the direct 2D model, reads the resulting simulated temperature file, and aligns it with observed data.
        - The simulation output file is renamed and saved in the appropriate output directory.
        - Assumes existence of helper functions: `modify_zone_parameters`, `run_direct_model_2D`.
        - Requires `os`, `pandas as pd`, and other relevant imports.
    """
    # Select the row corresponding to id_sim
    # Calculer le MSE pour chaque capteur sur les 8 premières valeurs et stocker dans une DataFrame
    row = param_table.loc[param_table['index_sim'] == id_sim]
    param_table_simul = row.copy()
    value_parameter = row.iloc[0].copy()
    # Modify zone parameters using the dictionary
    modify_zone_parameters(param_struct, zones_to_invert, parameters_to_invert, value_parameter)
    run_direct_model_2D(dir_ginette)
    # add "Time" tp name_sensor vector
    col_name= ["Time"] + name_sensor
    sim_temp = pd.read_csv("S_temp_PT100_t.dat", sep='\s+', header=0,names=col_name)
      # Ajouter une colonne de dates en utilisant la date de début de simulation
    date_simul_bg = pd.to_datetime(date_simul_bg)
    sim_temp['dates'] = date_simul_bg + pd.to_timedelta(sim_temp['Time'], unit='s')
    path_output_tempsim = os.path.join(repertory, f"OUTPUT_{Station}")
    new_file_name = f"{id_sim}_S_temp_PT100_t.dat"
    # save in the folder OUTPUT path_output_tempsim
    new_file_path = os.path.join(path_output_tempsim, new_file_name)
    os.rename("S_temp_PT100_t.dat", new_file_path)


    return param_table_simul
#_______________________________________________________________________________________
#
def analysis_gridsearch_2D(Station,repertory,date_simul_bg,obs_temp,name_sensor,param_table,id_sim,param_struct,zones_to_invert,parameters_to_invert,param_table_simul,sigma):
    """
    Calculate the mse for each sensor and store the results in a DataFrame.

    Parameters:
        Station (str): Name of the station.
        repertory (str): Directory where the simulation files are located.
        date_simul_bg (str): Start date of the simulation in 'YYYY-MM-DD' format.
        obs_temp (pd.DataFrame): DataFrame containing observed temperature data.
        name_sensor (list): List of sensor names to be used in the simulation.
        param_table (pd.DataFrame): DataFrame containing parameter values for the simulation.
        id_sim (int): Simulation ID.
        param_struct (list): List of parameter names in the structure.
        zones_to_invert (list): List of zones to be inverted.
        parameters_to_invert (list): List of parameters to be inverted.
        param_table_simul (pd.DataFrame): DataFrame to store simulation results.
        sigma (float): Standard deviation for the simulation.

    Returns:
        mse_df (pd.DataFrame): DataFrame containing MSE values for each sensor.
        param_table_simul (pd.DataFrame): Updated DataFrame with simulation parameters.

    Notes:
        - Modifies zone parameters using the provided structure and inversion lists.
        - Reads the resulting simulated temperature file and aligns it with observed data.
        - Removes the first two days of data from both simulated and observed datasets before MSE calculation.
        - Assumes existence of helper functions: `modify_zone_parameters`, `run_direct_model_2D`, `remove_first_two_days`, and `mse`.
        - Requires `os`, `pandas as pd`, and other relevant imports.
    """
    # Initialize dictionary to store MSE values for each sensor
    mse_values = {}
    
    # Create empty DataFrame to store MSE results with sensor names as columns
    # This will hold the calculated MSE for each temperature sensor
    mse_df = pd.DataFrame(columns=name_sensor)
    
    # STEP 1: Extract simulation parameters for the current simulation ID
    # Find the row in param_table that corresponds to our simulation ID
    row = param_table.loc[param_table['index_sim'] == id_sim]
    
    # Create a copy of the parameters for this simulation
    param_table_simul = row.copy()
    
    # Extract parameter values as a Series for easy access
    value_parameter = row.iloc[0].copy()
    
    # STEP 2: Define column names for the simulated temperature file
    # The simulation output file has Time column + 8 temperature sensors (Temp_1 to Temp_8)
    col_name = ["Time"] + [f"Temp_{i}" for i in range(1, 9)]
    
    # STEP 3: Construct the path to the simulation output file
    # Build the path to the OUTPUT folder for this station
    path_output_tempsim = os.path.join(repertory, f"OUTPUT_{Station}")
    
    # The simulation file is named with the simulation ID as prefix
    file_name = f"{id_sim}_S_temp_PT100_t.dat"
    
    # Complete file path to the simulation results
    file_path = os.path.join(path_output_tempsim, file_name)
    
    # STEP 4: Read the simulated temperature data
    # Load the simulation results file with space-separated values
    sim_temp = pd.read_csv(file_path, sep='\s+', header=0, names=col_name)
    
    # STEP 5: Add datetime column to simulation data
    # Convert the start date string to datetime object
    date_simul_bg = pd.to_datetime(date_simul_bg)
    
    # Create a 'dates' column by adding simulation time (in seconds) to start date
    # This allows us to align simulated and observed data by timestamp
    sim_temp['dates'] = date_simul_bg + pd.to_timedelta(sim_temp['Time'], unit='s')
    
    # STEP 6: Prepare observed data for comparison
    # Reset the index of observed temperature data to ensure clean indexing
    obs_temp.reset_index(drop=True, inplace=True)
    
    # Remove the first 2 days from simulation data to allow model stabilization
    # This is common practice to avoid initialization artifacts in numerical simulations
    sim_temp = remove_first_days_sim(sim_temp, 2)

    # Check size  of simulated data and obs data
    # Check if both datasets have the same number of rows
    if len(sim_temp) != len(obs_temp):
        print(f"Warning: Size mismatch - Simulated data: {len(sim_temp)} rows, Observed data: {len(obs_temp)} rows")
        
        # Find common dates between simulated and observed data
        # Convert observed data index to datetime if it's not already
        if not isinstance(obs_temp.index, pd.DatetimeIndex):
            obs_temp_dates = pd.to_datetime(obs_temp.index)
        else:
            obs_temp_dates = obs_temp.index
        
        # Find intersection of dates
        common_dates = obs_temp_dates.intersection(sim_temp['dates'])
        print(f"Found {len(common_dates)} common dates")
        
        if len(common_dates) > 0:
            # Filter both datasets to common dates
            sim_temp = sim_temp[sim_temp['dates'].isin(common_dates)]
            obs_temp = obs_temp[obs_temp_dates.isin(common_dates)]
            print(f"Both datasets filtered to {len(common_dates)} common dates")
        else:
            print("No common dates found, trimming to minimum length")
            min_length = min(len(sim_temp), len(obs_temp))
            sim_temp = sim_temp.iloc[:min_length]
            obs_temp = obs_temp.iloc[:min_length]
            print(f"Both datasets trimmed to {min_length} rows for consistency")
    else:
        print(f"Data sizes match: {len(sim_temp)} rows in both datasets")




    # STEP 7: Calculate MSE for each sensor
    # Loop through each temperature sensor to compute Mean Squared Error
    for col in name_sensor:  # Example: name_sensor = ['Temp1', 'Temp2', 'Temp3']
        
        # Extract simulated temperatures for this sensor
        prediction = sim_temp[col]
        
        # Extract observed temperatures for this sensor
        observation = obs_temp[col]
        
        # Calculate MSE between simulated and observed temperatures
        # The sigma parameter is used for normalization in the MSE calculation
        mse_val = mse(prediction, observation, sigma=sigma)
        
        # Store the MSE value in the results DataFrame
        mse_df.at[0, col] = mse_val
        
        # Optional: Print MSE for debugging (currently commented out)
        # print(f"Sensor {col}: MSE = {mse_val}")

    # STEP 8: Calculate total MSE across all sensors
    # Sum MSE values across all sensors to get overall model performance
    # This gives us a single metric to evaluate how well the model fits all observations
    mse_df['Total_mse'] = mse_df.sum(axis=1)

    # STEP 9: Return results
    # Return both the MSE DataFrame and the simulation parameters
    # mse_df contains MSE for each sensor plus total MSE
    # param_table_simul contains the parameters used for this simulation
    return mse_df, param_table_simul


#______________________________________________________________________________________
# function plot_marginal_1D

def plot_marginal_1D(df, param_name, param_units=None):
    """
    Plot the 1D marginal distribution P(param) for a given parameter, with optional unit in the legend.

    Args:
        df (pd.DataFrame): DataFrame containing parameter columns and a 'posterior' column.
        param_name (str): Name of the parameter column to plot (e.g., "k4").
        param_units (dict, optional): Dictionary mapping parameter names to their units (e.g., {"k4": "m/s"}).
    """
    grouped = df.groupby(param_name)["posterior"].sum()

    plt.figure(figsize=(6, 4))
    plt.plot(grouped.index, grouped.values, marker='o', label=param_name)

    # Prepare label with unit if provided
    if param_units and param_name in param_units:
        unit = param_units[param_name]
        legend_label = f"{param_name} [{unit}]"
        plt.xlabel(legend_label)
        plt.legend([legend_label])
    else:
        plt.xlabel(param_name)
        plt.legend([param_name])

    plt.ylabel(f"P({param_name})")
    plt.title(f"Marginal 1D: {param_name}")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

#_______________________
def plot_joint_posterior(df, param_cols, filename="grid_mse.png", param_units=None, cmap="viridis_r", sigma=None):
    """
    Plot 1D and 2D marginal MSE distributions from a grid search optimized for parameter sensitivity analysis.
    Specifically designed for permeability, porosity, and thermal conductivity sensitivity analysis.
    
    Args:
        df (DataFrame): DataFrame containing param_cols + 'Total_mse' columns
        param_cols (list): List of parameter column names to visualize
        filename (str): Output filename
        param_units (dict): Parameter units (e.g., {'k': 'm²', 'n': '-', 'l': 'W/m/K'})
        cmap (str): Colormap for heatmaps (viridis_r so low MSE appears yellow)
        sigma (float): Temperature uncertainty used in MSE calculation (for interpretation)
    
    Returns:
        dict: Dictionary containing best model parameters, sensitivity metrics, and statistics
    """
    # Initialize results dictionary
    results = {}

    # Find best model (minimum MSE)
    best_idx = df['Total_mse'].idxmin()
    best_model = df.loc[best_idx]
    results['best_model'] = best_model.to_dict()
    results['best_mse'] = best_model['Total_mse']

    # Calculate variance and standard deviation of Total_mse
    total_mse_values = df['Total_mse']
    mse_mean = total_mse_values.mean()
    mse_variance = total_mse_values.var()
    mse_std = total_mse_values.std()
    mse_sum= total_mse_values.sum()

    results['mse_statistics'] = {
        'mean': mse_mean,
        'variance': mse_variance,
        'std': mse_std,
        'sum' : mse_sum,
        'min': total_mse_values.min(),
        'max': total_mse_values.max()
    }

    # Calculate confidence intervals
    confidence_1sigma = mse_mean + mse_std
    confidence_2sigma = mse_mean + 2*mse_std
    models_within_1sigma = df[df['Total_mse'] <= confidence_1sigma]
    models_within_2sigma = df[df['Total_mse'] <= confidence_2sigma]
    # Seuils d'acceptabilité (modèles dans les X% meilleurs)
    threshold_10_percent = np.percentile(total_mse_values, 10)
    threshold_5_percent = np.percentile(total_mse_values, 5)
    models_within_10p = df[df['Total_mse'] <= threshold_10_percent]
    models_within_5p = df[df['Total_mse'] <= threshold_5_percent]

    print(f"Best model MSE: {results['best_mse']:.4f}")
    print(f"MSE variance: {mse_variance:.4f}")
    print(f"1σ confidence level: {confidence_1sigma:.4f}")
    print(f"2σ confidence level: {confidence_2sigma:.4f}")
    print(f"Seuil 10% meilleurs: {threshold_10_percent:.0f}")
    print(f"Seuil 5% meilleurs: {threshold_5_percent:.0f}")
    print(f"Models within 5p: {len(models_within_5p)}")


    # Create subplot grid
    n = len(param_cols)
    fig, axes = plt.subplots(n, n, figsize=(4*n, 4*n), squeeze=False)

    # Add colorbar axis
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar_drawn = False

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            param_x = param_cols[j]
            param_y = param_cols[i]
            
            # Prepare axis labels with units if provided
            xlabel = param_x
            ylabel = param_y
            if param_units:
                if param_x in param_units:
                    xlabel = f"{param_x} [{param_units[param_x]}]"
                if param_y in param_units:
                    ylabel = f"{param_y} [{param_units[param_y]}]"
            
            if i == j:  # Diagonal: 1D marginal plots
                # Group by parameter and sum MSE (inverse for probability-like visualization)
             # MODIFICATION : Grouper par paramètre et calculer statistiques Total_mse
                grouped_stats = df.groupby(param_x)['Total_mse'].agg(['sum', 'mean', 'min', 'std', 'count'])


                # Plot MSE curve
                ax.plot(grouped_stats.index, grouped_stats['sum'], 'b-', linewidth=2, label='MSE (sum)')
                # ERREUR : Vous ajoutez 5% au minimum au lieu de prendre 5% de plus que le minimum
                # Plage symétrique : minimum ± X%
                margin_5_percent = 0.05 * grouped_stats['sum'].min()
                margin_1_percent = 0.01 * grouped_stats['sum'].min()

                within_5_percent_lower = grouped_stats['sum'].min() - margin_5_percent
                within_5_percent_upper = grouped_stats['sum'].min() + margin_5_percent

                within_1_percent_lower = grouped_stats['sum'].min() - margin_1_percent
                within_1_percent_upper = grouped_stats['sum'].min() + margin_1_percent

                # Indices dans les plages symétriques
                idx_5 = (grouped_stats['sum'] >= within_5_percent_lower) & (grouped_stats['sum'] <= within_5_percent_upper)
                idx_1 = (grouped_stats['sum'] >= within_1_percent_lower) & (grouped_stats['sum'] <= within_1_percent_upper)
                                # Plot the acceptable region
                acceptable_params = grouped_stats.index[idx_5]
                acceptable_mse = grouped_stats['sum'][idx_5]
                ax.fill_between(acceptable_params, acceptable_mse, alpha=0.3, color='orange', 
                               label='Within 5% of best')
                # Plot the acceptable region
                acceptable_params = grouped_stats.index[idx_1]
                acceptable_mse = grouped_stats['sum'][idx_1]
                ax.fill_between(acceptable_params, acceptable_mse, alpha=0.3, color='green', 
                               label='Within 1% of best')


                # Mark best model
                best_param_val = best_model[param_x]
                best_mse_val = grouped_stats.loc[best_param_val, 'sum']
                ax.plot(best_param_val, best_mse_val, 'ro', markersize=8, label='Best model')
                

                ax.set_xlabel(xlabel)
                ax.set_ylabel('MSE')
                ax.grid(True, alpha=0.3)
                ax.legend(fontsize=8)
                
            else:  # Off-diagonal: 2D marginal plots
                # Create pivot table for heatmap
                pivot = df.pivot_table(index=param_y, columns=param_x, values='Total_mse', aggfunc='sum')
                
                # Create heatmap
                if not cbar_drawn:
                    im = sns.heatmap(pivot, ax=ax, cmap=cmap, cbar=True, cbar_ax=cbar_ax, 
                                   cbar_kws={'label': 'MSE'})
                    cbar_drawn = True
                else:
                    sns.heatmap(pivot, ax=ax, cmap=cmap, cbar=False)
                
                # Mark best model point
                best_x_val = best_model[param_x]
                best_y_val = best_model[param_y]
                
                # Find position in heatmap coordinates
                x_pos = list(pivot.columns).index(best_x_val) + 0.5
                y_pos = list(pivot.index).index(best_y_val) + 0.5
                
                ax.plot(x_pos, y_pos, 'ro', markersize=10, markeredgecolor='white', 
                       markeredgewidth=2, label='Best model')
                
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)



    plt.tight_layout(rect=[0, 0, 0.9, 0.92])
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.show()

    # Print detailed best model parameters
    print("\n" + "="*50)
    print("BEST MODEL PARAMETERS:")
    print("="*50)
    for param in param_cols:
        value = best_model[param]
        unit = param_units.get(param, '') if param_units else ''
        unit_str = f" [{unit}]" if unit else ""
        print(f"{param}: {value:.4f}{unit_str}")

    results['sensitivity_analysis'] = {}
    for param in param_cols:
        param_values = df[param].unique()
        param_mse_range = df.groupby(param)['Total_mse'].min().max() - df.groupby(param)['Total_mse'].min().min()
        results['sensitivity_analysis'][param] = {
            'range': param_mse_range,
            'normalized_sensitivity': param_mse_range / mse_std
        }

    print("\n" + "="*50)
    print("PARAMETER SENSITIVITY ANALYSIS:")
    print("="*50)
    for param, metrics in results['sensitivity_analysis'].items():
        print(f"{param}: MSE range = {metrics['range']:.4f}, Normalized sensitivity = {metrics['normalized_sensitivity']:.2f}")
    return results


# Plot joint posterior by sensor
def plot_joint_posterior_by_sensor(df, param_cols, col, param_units=None, cmap="viridis", filename="ind_posterior.png", plot_title=None):

    n = len(param_cols)
    # Find best model (minimum MSE)
    best_idx = df[col].idxmin()
    best_model = df.loc[best_idx]
    fig, axes = plt.subplots(n, n, figsize=(3*n, 3*n), squeeze=False)
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  
    cbar_drawn = False

    for i in range(n):
        for j in range(n):
            param_x = param_cols[j]
            param_y = param_cols[i]
            ax = axes[i, j]
            xlabel = param_cols[j]
            ylabel = param_cols[i]

            if param_units:
                if param_cols[j] in param_units:
                    xlabel = f"{param_cols[j]} [{param_units[param_cols[j]]}]"
                if param_cols[i] in param_units:
                    ylabel = f"{param_cols[i]} [{param_units[param_cols[i]]}]"

            if i == j:
                grouped = df.groupby(param_cols[i])[col].sum()
                                # Plage symétrique : minimum ± X%
                margin_5_percent = 0.05 * grouped.min()
                margin_1_percent = 0.01 * grouped.min()

                within_5_percent_lower = grouped.min() - margin_5_percent
                within_5_percent_upper = grouped.min() + margin_5_percent

                within_1_percent_lower = grouped.min() - margin_1_percent
                within_1_percent_upper = grouped.min() + margin_1_percent

                # Indices dans les plages symétriques
                idx_5 = (grouped >= within_5_percent_lower) & (grouped <= within_5_percent_upper)
                idx_1 = (grouped >= within_1_percent_lower) & (grouped <= within_1_percent_upper)
                                # Plot the acceptable region
                acceptable_params = grouped.index[idx_5]
                acceptable_mse = grouped[idx_5]
                ax.fill_between(acceptable_params, acceptable_mse, alpha=0.3, color='orange', 
                               label='Within 5% of best')
                # Plot the acceptable region
                acceptable_params = grouped.index[idx_1]
                acceptable_mse = grouped[idx_1]
                ax.fill_between(acceptable_params, acceptable_mse, alpha=0.3, color='green', 
                               label='Within 1% of best')
                # Mark best model
                best_param_val = best_model[param_x]
                best_mse_val = grouped.loc[best_param_val]
                ax.plot(best_param_val, best_mse_val, 'ro', markersize=8, label='Best model')
                
                ax.plot(grouped.index, grouped.values, color="black")
                ax.set_ylabel(f"MSE (sum) ({ylabel})")
                ax.set_xlabel(xlabel)
            else:
                pivot = df.pivot_table(index=param_cols[i], columns=param_cols[j], values=col, aggfunc='sum')
                if not cbar_drawn:
                    sns.heatmap(pivot, ax=ax, cmap=cmap, cbar=True, cbar_ax=cbar_ax)
                    cbar_drawn = True
                else:
                    sns.heatmap(pivot, ax=ax, cmap=cmap, cbar=False)
                # Mark best model point
                best_x_val = best_model[param_x]
                best_y_val = best_model[param_y]
                
                # Find position in heatmap coordinates
                x_pos = list(pivot.columns).index(best_x_val) + 0.5
                y_pos = list(pivot.index).index(best_y_val) + 0.5
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)

    if plot_title:
        plt.subplots_adjust(top=0.92)  # ajuster l'espace pour le titre

    if cbar_drawn:
        cbar_ax.set_ylabel("MSE (sum)", rotation=270, labelpad=15) 
    plt.tight_layout(rect=[0, 0, 0.9, 1])
    plt.savefig(filename)
    plt.show()
