# Direct_model.py
''' This module contains functions for setting up and running the Ginette model, 
as well as processing simulation results. The functions handle various tasks 
such as setting up model parameters, initial conditions, boundary conditions, 
and generating zone parameters. Additionally, the module includes utilities 
for running the simulation and processing its output.
Functions:
- setup_ginette_perm: Sets up the Ginette model parameters for steady-state simulations.
- setup_ginette_perm_2D: Configures the Ginette model for 2D steady-state simulations.
- setup_ginette: Sets up the Ginette model parameters for transient simulations.
- initial_conditions_perm_2D: Configures initial conditions for 2D steady-state simulations.
- initial_conditions: Sets up initial conditions for the Ginette model using input data.
- format_value: Formats a floating-point value in the format 00000000d+00.
- boundary_conditions_2D: Applies boundary conditions for 2D simulations and updates parameter files.
- boundary_conditions: Applies boundary conditions for transient simulations and updates parameter files.
- boundary_conditions_perm: Sets boundary conditions for permeability simulations.
- generate_zone_parameters: Writes zone parameter files based on the number of zones and parameter values.
- generate_zone_parameters_undef: Writes zone parameter files for undefined configurations.
- run_direct_model: Runs the Ginette model simulation and processes the output temperature data.
- remove_first_two_days: Removes the first four days from simulation and observation data.
- remove_first_days_sim: Removes the first few days from simulation data.
- reuse_end_in_inital: Copies the last simulation state to the initial condition file.
- The module assumes the presence of specific input files in the working directory.
- The functions rely on pandas and numpy for data manipulation and processing.
- The Ginette model executable must be available in the working directory for simulations. '''
import pandas as pd
import numpy as np
import subprocess
import os

def format_value(value):
    """
    Formate une valeur flottante dans le format 00000000d+00.
    """
    return "{:0=+12.2e}".format(value).replace('e', 'd')

def setup_ginette_perm(dt, state, nb_day, z_top, z_bottom, az, dz, date_simul_bg, dz_obs):
    """
    Sets up the Ginette model parameters and writes them to the appropriate files in steady state.
    Parameters:
    dt (float): Time step for the simulation.
    state (int): State of the simulation.
    nb_day (float): Number of days for the simulation.
    z_top (float): Top boundary of the model domain.
    z_bottom (float): Bottom boundary of the model domain.
    az (float): Total height of the model domain.
    dz (float): Cell height in the model domain.
    date_simul_bg (str): Start date of the simulation.
    dz_obs (float): Observation depth interval.
    Returns:
    list: A list of observation depths.
    """
    
    print("la simulation commence à", date_simul_bg)
    # number of cell
    # number of cell
    nb_cell=az/dz
    #-----------------------------------------------------------------
    ## write the setup of the modeled domain
    f_param_bck = open("E_parametre_bck.dat", "r")
    f_param_new = open("E_parametre.dat", 'w')
    setup_model = f_param_bck.read()
    setup_model = setup_model.replace('[dt]', '%06.0fD+00' % dt)
    setup_model = setup_model.replace('[state]', '%1i' % state)
    setup_model = setup_model.replace('[nb_day]', '%06.0f' % nb_day)
    setup_model = setup_model.replace('[z_top]', '%7.3e' % z_top)
    setup_model = setup_model.replace('[z_bottom]', '%7.2e' % z_bottom)
    setup_model = setup_model.replace('[az]', '%7.3e' % az)
    setup_model = setup_model.replace('[dz]', '%6.2e' % dz)
    setup_model = setup_model.replace('[nb_cell]', '%05.0f' % nb_cell)
    setup_model= setup_model.replace('[itsortie]', '%08.0f' % dt)

    # Observation positions x 0.50000
    # Observation in meter
    Obs1 = z_top - dz_obs
    Obs2 = z_top - dz_obs * 2
    Obs3 = z_top - dz_obs * 3
    Obs4 = z_top - dz_obs * 4

    # Create z_obs vector
    z_obs = [Obs1, Obs2, Obs3, Obs4]
    ## write the parameters
    cell1 = abs(Obs1 / dz)
    cell2 = abs(Obs2 / dz)
    cell3 = abs(Obs3 / dz)
    cell4 = abs(Obs4 / dz)
    setup_model = setup_model.replace('[cell1]', '%05.0f' % cell1)
    setup_model = setup_model.replace('[cell2]', '%05.0f' % cell2)
    setup_model = setup_model.replace('[cell3]', '%05.0f' % cell3)
    setup_model = setup_model.replace('[cell4]', '%05.0f' % cell4)
    f_param_new.write(setup_model)
    f_param_bck.close()
    f_param_new.close()
    
    #-----------------------------------------------------------------
    f_param_therm=open("E_p_therm_bck.dat", "r")
    f_param_therm_new=open("E_p_therm.dat", "w")
    therm=f_param_therm.read()
    therm = therm.replace('[state]', '%1i' % state)
    f_param_therm_new.write(therm)
    f_param_therm_new.close()
    f_param_therm.close()
    

    return z_obs
    
    #-----------------------------------------------------------------


def setup_ginette_perm_2D(pt100_coord,nb_cell,nb_col,nb_row):
    """
    Sets up the Ginette model parameters and writes them to the appropriate files in steady state.
    Parameters:
    dt (float): Time step for the simulation.
    state (int): State of the simulation.
    nb_day (float): Number of days for the simulation.
    date_simul_bg (str): Start date of the simulation.
    pt100_coord (DataFrame): DataFrame containing the coordinates of the PT100 sensors.
                            hobo  pt100    x       z  index  distance  xmaille  zmaille  n_maille
    Returns:
    """
    

    # time step
    dt=900
    # state
    state=0
    nb_day=10
    itsortie=dt 
    #-----------------------------------------------------------------
    ## write the setup of the modeled domain
    f_param_bck = open("E_parametre_bck.dat", "r")
    f_param_new = open("E_parametre.dat", 'w')
    setup_model = f_param_bck.read()
    setup_model = setup_model.replace('[dt]', '%06.0fD+00' % dt)
    setup_model = setup_model.replace('[nb_day]', '%06.0f' % nb_day)    
    setup_model = setup_model.replace('[state]', '%1i' % state)
    setup_model = setup_model.replace('[nb_cell]', '%05.0f' % nb_cell)
    setup_model = setup_model.replace('[nb_col]', '%05.0f' % nb_col)
    setup_model = setup_model.replace('[nb_row]', '%05.0f' % nb_row)
    setup_model= setup_model.replace('[itsortie]', '%08.0f' % dt)

    # read pt100_coord to get the number of rows 
    nb_pt100=len(pt100_coord)
    print(pt100_coord.head())
    # number cell to replace = nb_pt100
    # create the int of cell1, cell2, cell3, cell4,cell5, cell6,cell7, cell8 =10
    for i in range(8):
        # create the name of variable cell1, cell2, cell3, cell4,cell5, cell6,cell7, cell8 in function of i
        cell_name = 'cell' + str(i+1)
        # get the value of cell1, cell2, cell3, cell4,cell5, cell6,cell7, cell8 in function of i
        if i < nb_pt100:
            cell_value = pt100_coord.iloc[i]['n_maille']
        else:
            cell_value = 10
        setup_model = setup_model.replace('[' + cell_name + ']', '%05.0f' % cell_value)
        # write the value of cell1, cell2, cell3, cell4,cell5, cell6,cell7, cell8 in the file
        setup_model = setup_model.replace('[' + cell_name + ']', '%05.0f' % cell_value)



    f_param_new.write(setup_model)
    f_param_bck.close()
    f_param_new.close()
    
    #-----------------------------------------------------------------
    f_param_therm=open("E_p_therm_bck.dat", "r")
    f_param_therm_new=open("E_p_therm.dat", "w")
    therm=f_param_therm.read()
    therm = therm.replace('[state]', '%1i' % state)
    f_param_therm_new.write(therm)
    f_param_therm_new.close()
    f_param_therm.close()
    

     


def setup_ginette(dt, state, nb_day, z_top, z_bottom, az, dz, date_simul_bg,dz_obs):
    """
    Sets up the Ginette model parameters and writes them to the appropriate files in transient state.
    Parameters:
    dt (float): Time step for the simulation.
    state (int): State of the simulation.
    nb_day (float): Number of days for the simulation.
    z_top (float): Top boundary of the model domain.
    z_bottom (float): Bottom boundary of the model domain.
    az (float): Total height of the model domain.
    dz (float): Cell height in the model domain.
    date_simul_bg (str): Start date of the simulation.
    dz_obs (float): Observation depth interval.
    Returns:
    list: A list of observation depths.
    """
    print("la simulation commence à", date_simul_bg)
    
    # number of cell
    nb_cell=az/dz
    #-----------------------------------------------------------------
    ## write the setup of the modeled domain
    f_param_bck = open("E_parametre_bck.dat", "r")
    f_param_new = open("E_parametre.dat", 'w')
    setup_model = f_param_bck.read()
    setup_model = setup_model.replace('[dt]', '%06.0fD+00' % dt)
    setup_model = setup_model.replace('[state]', '%1i' % state)
    setup_model = setup_model.replace('[nb_day]', '%06.0f' % nb_day)
    setup_model = setup_model.replace('[z_top]', '%7.3e' % z_top)
    setup_model = setup_model.replace('[z_bottom]', '%7.2e' % z_bottom)
    setup_model = setup_model.replace('[az]', '%7.3e' % az)
    setup_model = setup_model.replace('[dz]', '%6.2e' % dz)
    setup_model = setup_model.replace('[nb_cell]', '%05.0f' % nb_cell)
    setup_model= setup_model.replace('[itsortie]', '%08.0f' % dt)

    # Observation positions x 0.50000
    # Observation in meter
    Obs1 = z_top - dz_obs
    Obs2 = z_top - dz_obs * 2
    Obs3 = z_top - dz_obs * 3
    Obs4 = z_top - dz_obs * 4

    # Create z_obs vector
    z_obs = [Obs1, Obs2, Obs3, Obs4]
    ## write the parameters
    cell1 = abs(Obs1 / dz)
    cell2 = abs(Obs2 / dz)
    cell3 = abs(Obs3 / dz)
    cell4 = abs(Obs4 / dz)
    setup_model = setup_model.replace('[cell1]', '%05.0f' % cell1)
    setup_model = setup_model.replace('[cell2]', '%05.0f' % cell2)
    setup_model = setup_model.replace('[cell3]', '%05.0f' % cell3)
    setup_model = setup_model.replace('[cell4]', '%05.0f' % cell4)
    f_param_new.write(setup_model)
    f_param_bck.close()
    f_param_new.close()
    
    #-----------------------------------------------------------------
    f_param_therm=open("E_p_therm_bck.dat", "r")
    f_param_therm_new=open("E_p_therm.dat", "w")
    therm=f_param_therm.read()
    therm = therm.replace('[state]', '%1i' % state)
    f_param_therm_new.write(therm)
    f_param_therm_new.close()
    f_param_therm.close()
    
    #-----------------------------------------------------------------

    return z_obs
    
    #-----------------------------------------------------------------
def initial_conditions_perm_2D():
    """
    Sets up the initial conditions for a 2D model by reading a template file,
    replacing placeholders with specific values, and writing the modified
    content to a new file.
    This function performs the following steps:
    1. Opens a template file named "E_cdt_initiale_bck.dat" in read mode.
    2. Reads the content of the template file.
    3. Replaces placeholders '[chg_i]' and '[temp_i]' in the content with
       specific integer values (0 in this case).
    4. Writes the modified content to a new file named "E_cdt_initiale.dat".
    5. Closes both the input and output files.
    Note:
    - The placeholders '[chg_i]' and '[temp_i]' in the template file must
      exist for the replacement to work correctly.
    - The function assumes that the files are located in the current working
      directory.
    Raises:
    - FileNotFoundError: If the template file "E_cdt_initiale_bck.dat" does
      not exist.
    - IOError: If there are issues reading from or writing to the files.
    """

    ## write the setup of the modeled domain
    f_IC_bck = open("E_cdt_initiale_bck.dat", "r")
    f_IC_new = open("E_cdt_initiale.dat", 'w')
    setup_IC=f_IC_bck.read()
    ichg_i=0
    itemp_i=0
    setup_IC = setup_IC.replace('[chg_i]', '%1i' % ichg_i)
    setup_IC = setup_IC.replace('[temp_i]', '%1i' % itemp_i)
    f_IC_new.write(setup_IC)
    f_IC_bck.close()
    f_IC_new.close()   


def initial_conditions(all_data, z_top, z_bottom, dz, z_obs):
    """
    Sets up the initial conditions for the Ginette model.
    Parameters:
    all_data (DataFrame): DataFrame containing the initial data.
    z_top (float): Top boundary of the model domain.
    z_bottom (float): Bottom boundary of the model domain.
    dz (float): Cell height in the model domain.
    z_obs (list): List of observation depths.
    """
    # Initial conditions
    with open("E_temperature_initiale.dat", "w") as f_temp_IC:
        # Dynamically select temperature columns based on the names in all_data
        if 'T_top' in all_data.columns and 'T_bottom' in all_data.columns:
            temp_columns = ['T_top', 'T_bottom']
            z_temps = np.array([z_top + 0.005, z_bottom])
        elif 'TempMolo' in all_data.columns and 'Temp4' in all_data.columns:
            temp_columns = ['TempMolo'] + [f'Temp{i+1}' for i in range(len(z_obs))]
            z_temps = np.array([z_top + 0.005] + z_obs)
        else:
            print("Error: No initial temperature data found.")
            return  

        initial_temps = all_data.iloc[0][temp_columns]
        initial = pd.DataFrame({'z': z_temps, 'T': initial_temps})
        initial['z'] = initial['z'].astype(float)
        initial['T'] = initial['T'].astype(float)


        # Define the range of depths for interpolation
        z_values = np.arange(z_bottom + dz / 2, z_top + 0.001, dz)  # Depth range from z_bottom to z_top

        # Perform linear interpolation
        Temp_init = np.interp(z_values, initial['z'][::-1], initial['T'][::-1])  # Reversed to align with depths

        # Creating a DataFrame for the interpolated data
        interpolated_temp = pd.DataFrame({'z': z_values, 'T': Temp_init})
        interpolated_temp_sorted = interpolated_temp.sort_values(by='z', ascending=False)

        interpolated_temp_sorted['T'].to_csv(f_temp_IC, index=False, sep='\n', header=False)



  # Initial conditions for hydraulic head
    with open("E_charge_initiale.dat", "w") as f_chg_IC:
        if 'h_top' in all_data.columns and 'h_bottom' in all_data.columns:
            chg_columns = ['h_top', 'h_bottom']
            z_chg = np.array([z_top , z_bottom])
        elif 'deltaP' in all_data.columns:
            all_data['h_top'] = all_data['deltaP']
            all_data['h_bottom'] = 0
            chg_columns = ['h_top', 'h_bottom']
            z_chg = np.array([z_top, z_bottom])
        else:
            print("Error: No initial charge data found.")
            return

        # Perform linear interpolation
        initial_chg = all_data.iloc[0][chg_columns]
        initial = pd.DataFrame({'z': z_chg, 'chg': initial_chg})
        initial['z'] = initial['z'].astype(float)
        initial['chg'] = initial['chg'].astype(float)

        # Define the range of depths for interpolation
        z_values = np.arange(z_bottom + dz / 2, z_top, dz)  # Depth range from z_bottom to z_top

        # Perform linear interpolation
        chg_init = np.interp(z_values, initial['z'][::-1], initial['chg'][::-1])  # Reversed to align with depths

        # Creating a DataFrame for the interpolated data
        interpolated_chg = pd.DataFrame({'z': z_values, 'chg': chg_init})
        interpolated_chg_sorted = interpolated_chg.sort_values(by='z', ascending=False)

        interpolated_chg_sorted['chg'].to_csv(f_chg_IC, index=False, sep='\n', header=False)



def boundary_conditions_2D(all_data,dt):
    """
        Updates boundary condition parameters for a 2D model and writes them to files.
        This function reads a template file containing boundary condition parameters,
        modifies its content based on the provided data, and writes the updated
        parameters to a new file. Additionally, it saves specific columns of the
        input data to separate files for further processing.
        Args:
            all_data (pandas.DataFrame): A DataFrame containing the boundary condition
                data. Expected columns include:
                - 'H_RD': Hydraulic head at the right  boundary.
                - 'H_RG': Hydraulic head at the left  boundary.
                - 'T_RD': Temperature at the right boundary.
                - 'T_RG': Temperature at the left boundary.
                - 'T_RIV': Temperature at the river boundary.
                - 'H_RIV': Hydraulic head at the river boundary.
            dt (float): The current time step value, used to update the boundary
                condition file.
        File Outputs:
            - 'E_chargeT_RD.dat': Contains the 'H_RD' column values.
            - 'E_chargeT_RG.dat': Contains the 'H_RG' column values.
            - 'E_tempT_RD.dat': Contains the 'T_RD' column values.
            - 'E_tempT_RG.dat': Contains the 'T_RG' column values.
            - 'E_tempT_Riv.dat': Contains the 'T_RIV' column values.
            - 'E_chargeT_Riv.dat': Contains the 'H_RIV' column values.
            - 'E_cdt_aux_limites.dat': Updated boundary condition parameters.
        Notes:
            - The function reads from a template file named 'E_cdt_aux_limites_bck.dat'.
            - The template file must contain placeholders such as '[iclchgt]', '[itlecture]',
              '[valcl_haut]', and '[valcl_bas]', which are replaced with actual values.
            - The `format_value` function is assumed to be defined elsewhere in the code
              to format the 'top' and 'bot' values appropriately.
        Returns:
            None
    """
    f_param_lec=open("E_cdt_aux_limites_bck.dat", "r")
    f_param_lec_new=open("E_cdt_aux_limites.dat", "w")
    
    lec_bc=f_param_lec.read()
    # Save value of H_RD in 'E_chargeT_RD.dat'
    all_data[['H_RD']].to_csv('E_chargeT_RD.dat', sep=' ', index=False, header=False)

    # Save value of H_RG in 'E_chargeT_RG.dat'
    all_data[['H_RG']].to_csv('E_chargeT_RG.dat', sep=' ', index=False, header=False)

    # Save value of T_RD in 'E_tempT_RD.dat'
    all_data[['T_RD']].to_csv('E_tempT_RD.dat', sep=' ', index=False, header=False)

    # Save value of T_RG in 'E_tempT_RG.dat'
    all_data[['T_RG']].to_csv('E_tempT_RG.dat', sep=' ', index=False, header=False)

    # Save value of T_RIV in 'E_tempT_Riv.dat'
    all_data[['T_RIV']].to_csv('E_tempT_Riv.dat', sep=' ', index=False, header=False)

    # Save value of H_RIV in 'E_chargeT_Riv.dat'
    all_data[['H_RIV']].to_csv('E_chargeT_Riv.dat', sep=' ', index=False, header=False)

    iclchgt = 1
    lec_bc = lec_bc.replace('[iclchgt]', '%1i' % iclchgt)
    lec_bc = lec_bc.replace('[itlecture]', '%08.0f' % dt)
    lec_bc = lec_bc.replace('[valcl_haut]', format_value(all_data['H_RIV'].iloc[0]))
    lec_bc = lec_bc.replace('[valcl_gauche]', format_value(all_data['H_RG'].iloc[0]))
    lec_bc = lec_bc.replace('[valcl_droite]', format_value(all_data['H_RD'].iloc[0]))
    lec_bc = lec_bc.replace('[valclt_haut]', format_value(all_data['T_RIV'].iloc[0]))
    lec_bc = lec_bc.replace('[valclt_gauche]', format_value(all_data['T_RG'].iloc[0]))
    lec_bc = lec_bc.replace('[valclt_droite]', format_value(all_data['T_RD'].iloc[0]))
    lec_bc = lec_bc.replace('[itlecture]', '%08.0f' % dt) 
    
    # Écrire les modifications dans le fichier
    f_param_lec_new.write(lec_bc)
    return

def boundary_conditions(all_data,dt):

    """
    Apply boundary conditions to the given data and update parameter files in transient state.
    This function reads boundary condition parameters from a file, updates the 
    boundary conditions based on the provided data, and writes the updated 
    parameters back to a new file. It also saves the charge and temperature 
    boundary conditions to separate files.
    Parameters:
    all_data (DataFrame): A pandas DataFrame containing the boundary condition data.
                          Expected columns are 'deltaP', 'h_top', 'h_bottom', 
                              'T_top', 'T_bottom', 'TempMolo', and 'Temp4'.
    dt (float): A time step value used for updating the parameter file.
    Returns:
    str: The updated boundary condition parameters as a string.
    Raises:
    ValueError: If neither flow nor temperature boundary conditions are found in the data.
    Notes:
        - The function expects the input file "E_cdt_aux_limites_bck.dat" to be present in the 
          working directory.
        - The function writes the updated boundary condition parameters to "E_cdt_aux_limites.dat".
        - Charge boundary conditions are saved to "E_charge_t.dat".
        - Temperature boundary conditions are saved to "E_temp_t.dat".
    """
    
    f_param_lec=open("E_cdt_aux_limites_bck.dat", "r")
    f_param_lec_new=open("E_cdt_aux_limites.dat", "w")
    
    lec_bc=f_param_lec.read()
    
    if 'deltaP' in all_data.columns:    # Boundary conditions
        all_data['bot'] = 0
        all_data['top'] = all_data['deltaP']
    elif 'h_top' in all_data.columns and 'h_bottom' in all_data.columns:
        all_data['bot'] = all_data['h_bottom']
        all_data['top'] = all_data['h_top']
    else:
        print("Error: No flow boundary conditions found.")
        return
    
    # Save charge boundary conditions
    all_data[['top', 'bot']].to_csv('E_charge_t.dat', sep=' ', index=False, header=False)
    
    # Save temperature boundary conditions
    if 'T_top' in all_data.columns and 'T_bottom' in all_data.columns:
        all_data[['T_top', 'T_bottom']].to_csv('E_temp_t.dat', sep=' ', index=False, header=False)
        lec_bc = lec_bc.replace('[valclt_haut]', format_value(all_data['T_top'].iloc[0]))
        lec_bc = lec_bc.replace('[valclt_bas]', format_value(all_data['T_bottom'].iloc[0]))
    elif 'TempMolo' in all_data.columns and 'Temp4' in all_data.columns:
        all_data[['TempMolo', 'Temp4']].to_csv('E_temp_t.dat', sep=' ', index=False, header=False)
        lec_bc = lec_bc.replace('[valclt_haut]', format_value(all_data['TempMolo'].iloc[0]))
        lec_bc = lec_bc.replace('[valclt_bas]', format_value(all_data['Temp4'].iloc[0]))
    else:
        print("Error: No temperature boundary conditions found.")
        return


    iclchgt = 1
    lec_bc = lec_bc.replace('[iclchgt]', '%1i' % iclchgt)
    lec_bc = lec_bc.replace('[itlecture]', '%08.0f' % dt)
    lec_bc = lec_bc.replace('[valcl_haut]', format_value(all_data['top'].iloc[0]))
    lec_bc = lec_bc.replace('[valcl_bas]', format_value(all_data['bot'].iloc[0]))
    
    # Écrire les modifications dans le fichier
    f_param_lec_new.write(lec_bc)
    return


# 

#-----------------------------------------------------------------


def boundary_conditions_perm_2D(all_data):
    """
    Sets boundary conditions for permeability based on the provided data and writes them to a file.
    Parameters:
    all_data (DataFrame): A pandas DataFrame containing the necessary boundary condition data. 
                        Expected columns are 'deltaP', 'h_top', 'h_bottom', 'T_top', 'T_bottom', 
                        'TempMolo', and 'Temp4'.
    The function performs the following steps:
    1. Opens the existing boundary condition file "E_cdt_aux_limites_perm_bck.dat" for reading.
    2. Opens a new file "E_cdt_aux_limites.dat" for writing the updated boundary conditions.
    3. Reads the content of the existing boundary condition file.
    4. Checks for flow boundary conditions in the DataFrame:
    -H_RD, H_RIV,H_RG
    5. Checks for temperature boundary conditions in the DataFrame:
    -  'T_RIV' 
    -  'T_RG'
    -  'T_RD'
    - If neither condition is met, prints an error message and exits.
    6. Sets the charge boundary condition placeholder '[iclchgt]' to 0.
    7. Replaces placeholders '[valcl_haut]', '[valcl_gauche]' and '[valcl_droite]' in the file content with 'riv', 'RD' and 'RG' values.
    8. Writes the modified content to the new boundary condition file.
    9. Closes both the new and existing boundary condition files.
    Returns:
    None
    """   
        #-----------------------------------------------------------------
    f_param_lec=open("E_cdt_aux_limites_perm_bck.dat", "r")
    f_param_lec_new=open("E_cdt_aux_limites.dat", "w")
    lec_bc=f_param_lec.read()
    # Test if the three boundary conditions are present in the DataFrame
    if 'H_RIV' in all_data.columns and 'H_RD' in all_data.columns and 'H_RG' in all_data.columns:
        print('Flow boundaries columns names are OK')
    else:
        print("Error: Missing one or more flow boundary conditions (H_RIV, H_RD, H_RG).")
        return
    
    # Save charge boundary conditions

    
    # Save temperature boundary conditions
    if 'T_RIV' in all_data.columns and 'T_RG' in all_data.columns and 'T_RD' in all.data.columns:
        print('Heat boundaries columns names are Ok')
    else:
        print("Error: No temperature boundary conditions found.")
        return


    iclchgt = 0
    lec_bc = lec_bc.replace('[iclchgt]', '%1i' % iclchgt)
    lec_bc = lec_bc.replace('[valcl_haut]', '%70f' % all_data['H_RIV'].iloc[0])
    lec_bc = lec_bc.replace('[valcl_gauche]', '%70f' % all_data['H_RG'].iloc[0])
    lec_bc = lec_bc.replace('[valcl_droite]', '%70f' % all_data['H_RD'].iloc[0])
    lec_bc = lec_bc.replace('[valclt_haut]', '%70f' % all_data['T_RIV'].iloc[0])
    lec_bc = lec_bc.replace('[valclt_gauche]', '%70f' % all_data['T_RG'].iloc[0])
    lec_bc = lec_bc.replace('[valclt_droite]', '%70f' % all_data['T_RD'].iloc[0])
    # Écrire les modifications dans le fichier
    f_param_lec_new.write(lec_bc)

    # Fermer les fichiers
    f_param_lec_new.close()
    f_param_lec.close()

# 



   
def generate_zone_parameters(z_bottom, dz, nb_zone, alt_thk, REF_k, REF_n, REF_l, REF_r, REF_k2=None, REF_n2=None, REF_l2=None, REF_r2=None):
    """
    Writes the zone parameter file based on the number of zones and parameter values.

    Parameters:
    - z_bottom: Bottom boundary of the model domain.
    - dz: Cell height in the model domain.
    - nb_zone: Number of zones (1 or 2).
    - alt_thk: Altitude threshold for zone separation.
    - REF_k: Intrinsic permeability for zone 1.
    - REF_n: Porosity for zone 1.
    - REF_l: Thermal conductivity for zone 1.
    - REF_r: Density for zone 1.
    - REF_k2: Intrinsic permeability for zone 2 (if nb_zone == 2).
    - REF_n2: Porosity for zone 2 (if nb_zone == 2).
    - REF_l2: Thermal conductivity for zone 2 (if nb_zone == 2).
    - REF_r2: Density for zone 2 (if nb_zone == 2).
    """
    ########### Zone of parameters
    f_coor = open("E_coordonnee.dat", "w")
    f_zone = open("E_zone.dat", 'w')

    coord = pd.DataFrame()    

    # Coodrinate  
    zvalues = np.sort(np.arange(z_bottom + dz / 2, dz / 2, dz))
    xvalues = np.array([0.5])
    zz, xx = np.meshgrid(zvalues, xvalues)
    NT = np.prod(zz.shape)  # Remplacement de np.product par np.prod
    data = {
        "x": np.reshape(xx, NT),
        "z": np.reshape(zz, NT)
    }
    coord = pd.DataFrame(data=data)
    coord['id'] = coord.index.values.astype(int)
    coord['id'] = coord['id'] + 1
    cols = coord.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    coord = coord[cols] 
    coord.to_csv(f_coor, index=False, sep=' ', header=False)
    
    # zone parameter by cell (homogenous domain = 1 zone)
    coord['zone'] = 1


    # Pour plusieurs zones modification AR
    if nb_zone >= 2:
        coord['zone'] = np.where(coord['z'] >= alt_thk, 2, coord['zone'])
    
    # Write new ginette files
    coord.zone.to_csv(f_zone, index=False, header=False)

    # close files    
    f_zone.close()
    f_coor.close()
    
    param_zone = None
 
    if nb_zone == 1:
        k1= 10**REF_k
        # Lire le fichier de sauvegarde des paramètres de zone
        with open("E_zone_parameter_bck.dat", "r") as f_paramZ_bck:
            param_zone = f_paramZ_bck.read()
        # Remplacer les valeurs des paramètres
        param_zone = param_zone.replace('[k1]', '%8.2e' % k1)
        param_zone = param_zone.replace('[n1]', '%6.2f' % REF_n)
        param_zone = param_zone.replace('[l1]', '%6.2f' % REF_l)
        param_zone = param_zone.replace('[r1]', '%6.2f' % REF_r)
    elif nb_zone == 2:
        # Lire le fichier de sauvegarde des paramètres de zone pour 2 zones
        with open("E_zone_parameter_bck_2zones.dat", "r") as f_paramZ_bck:
            param_zone = f_paramZ_bck.read()
        k1= 10**REF_k
        k2= 10**REF_k2
        # Remplacer les valeurs des paramètres pour les deux zones
        param_zone = param_zone.replace('[k1]', '%8.2e' % k1)
        param_zone = param_zone.replace('[n1]', '%6.2f' % REF_n)
        param_zone = param_zone.replace('[l1]', '%6.2f' % REF_l)
        param_zone = param_zone.replace('[r1]', '%6.2f' % REF_r)
        param_zone = param_zone.replace('[k2]', '%8.2e' % k2)
        param_zone = param_zone.replace('[n2]', '%6.2f' % REF_n2)
        param_zone = param_zone.replace('[l2]', '%6.2f' % REF_l2)
        param_zone = param_zone.replace('[r2]', '%6.2f' % REF_r2)
    # Ouvrir le fichier de paramètres de zone
    f_paramZ_new = open("E_zone_parameter.dat", 'w')
    # Écrire les paramètres dans le nouveau fichier
    f_paramZ_new.write(param_zone)
    f_paramZ_new.close()

def zone_parameters_undef(nb_zone, parameters):
    """
    Generate a file "E_zone_parameter_bck.dat" with the parameters of the zones.

    Parameters:
    - nb_zone: Number of zones.
    - parameters: List of parameter names to be written in the file (e.g., ['k', 'n', 'l', 'cpm', 'r']).
    """
    with open("E_zone_parameter_bck.dat", "w") as f:
        for i in range(nb_zone):
            line = f"{i+1}"
            for param in parameters:
                line += f"\t[{param}{i+1}]"
            if i < nb_zone - 1:  # Avoid adding a newline at the end of the file
                f.write(line + "\n")
            else:
                f.write(line)


def generate_zone_parameters_undef(nb_zone,value_zone_parameter):
    """
    Writes the zone parameter file based on the number of zones and parameter values.

    Parameters:
    - nb_zone: 
    - value_zone_parameter: dataframe of parameter values for each zone, dataframe: with the 
    colum name are the same as parameters and each row is the value of the parameter for each zone.
    """
    f_param_zone=open("E_zone_parameter_bck.dat", "r")
    f_param_zone_new=open("E_zone_parameter.dat", "w")
    param_zone=f_param_zone.read()
    # if the parameter k exist in the file replace by k=10**k value_zone_parameter is a dataframe

    if 'k' in value_zone_parameter.columns:
        for i in range(nb_zone):
            # transform the value of k to 10**k
            value_zone_parameter = value_zone_parameter.copy()
            value_zone_parameter.loc[i, 'k'] = 10**value_zone_parameter.loc[i, 'k']


    # remplacer les valeurs des paramètres
    for i in range(nb_zone):
        for j in range(len(value_zone_parameter.columns)):
            param_zone = param_zone.replace('[' + value_zone_parameter.columns[j]
                                             + str(i+1) + ']', '{:0=12.2e}'.format
                                             (value_zone_parameter.iloc[i, j]).replace('e', 'd'))
    # write in new file E_zone_parameter.dat
    f_param_zone_new.write(param_zone)
    f_param_zone.close()
    f_param_zone_new.close()
    #-----------------------------------------------------------------          



def run_direct_model(date_simul_bg,z_bottom, dz, nb_zone, alt_thk, REF_k, REF_n, REF_l, REF_r, REF_k2=None, REF_n2=None, REF_l2=None, REF_r2=None):
    """
        Run the direct model simulation and process the output temperature data.
        Parameters:
        date_simul_bg (str): The start date of the simulation in a string format.
        z_bottom (float): The bottom depth of the simulation zone.
        dz (float): The depth increment for each zone.
        nb_zone (int): The number of zones in the simulation.
        alt_thk (float): The thickness of the alteration zone.
        REF_k (float): Reference parameter k for the simulation.
        REF_n (float): Reference parameter n for the simulation.
        REF_l (float): Reference parameter l for the simulation.
        REF_r (float): Reference parameter r for the simulation.
        REF_k2 (float, optional): Secondary reference parameter k for the simulation. Defaults to None.
        REF_n2 (float, optional): Secondary reference parameter n for the simulation. Defaults to None.
        REF_l2 (float, optional): Secondary reference parameter l for the simulation. Defaults to None.
        REF_r2 (float, optional): Secondary reference parameter r for the simulation. Defaults to None.
        Returns:
        pd.DataFrame: A DataFrame containing the simulation time, temperatures, and corresponding dates.
        Notes:
        - This function generates zone parameters and runs the 'ginette' simulation.
        - It reads the output temperature data from files and merges them into a single DataFrame.
        - The resulting DataFrame includes a 'dates' column calculated from the simulation start date and time.
    """


    generate_zone_parameters(z_bottom, dz, nb_zone, alt_thk, REF_k, REF_n, REF_l, REF_r, REF_k2, REF_n2, REF_l2, REF_r2)

    # run ginette
    subprocess.call(["./ginette"])
    # Colonnes pour les températures
    column_names = ['Temp_1', 'Temp_2', 'Temp_3']

    # Liste pour stocker les DataFrames
    all_dfs = []
    # Lire et stocker les DataFrames
    for i in range(1, 4):
        file_path = f'Sim_temperature_maille{i}_t.dat'
        df = pd.read_csv(file_path, sep='\s+', header=None, names=['Time', f'Temp{i}'])
        all_dfs.append(df)

    # Fusionner les DataFrames sur la colonne 'Time'
    sim_temp = all_dfs[0]  # Utiliser le premier DataFrame comme base
    for df in all_dfs[1:]:
        sim_temp = sim_temp.merge(df, on='Time')
    
      # Ajouter une colonne de dates en utilisant la date de début de simulation
    date_simul_bg = pd.to_datetime(date_simul_bg)
    sim_temp['dates'] = date_simul_bg + pd.to_timedelta(sim_temp['Time'], unit='s')

    return sim_temp

def run_direct_model_2D():
    """
        Run the direct model simulation and process the output temperature data.
        Parameters:
        date_simul_bg (str): The start date of the simulation in a string format.
        Returns:
        pd.DataFrame: A DataFrame containing the simulation time, temperatures, and corresponding dates.
        -S_temp_PT100_t.dat
        Notes:
        - This function generates zone parameters and runs the 'ginette' simulation.
        - It reads the output temperature data from files and merges them into a single DataFrame.
        - The resulting DataFrame includes a 'dates' column calculated from the simulation start date and time.
    """



    # run ginette
    subprocess.call(["./ginette"])
    # Colonnes pour les températures

    return 

def remove_first_two_days(sim_temp, obs_temp):
    """
    Removes the first four days from the sim_temp and obs_temp DataFrames.

    Parameters:
    - sim_temp: DataFrame containing simulated temperature data with a 'Time' column in seconds.
    - obs_temp: DataFrame containing observed temperature data with a datetime index.

    Returns:
    - sim_temp_filtered: Filtered sim_temp DataFrame.
    - obs_temp_filtered: Filtered obs_temp DataFrame.
    """
    # Filtrer sim_temp pour supprimer les deux premiers jours (86400 secondes par jour)
    sim_temp_filtered = sim_temp[sim_temp['Time'] >= 86400 * 4  ]

    # Convertir l'index de obs_temp en format datetime si ce n'est pas déjà fait
    obs_temp.index = pd.to_datetime(obs_temp.index)

    # Définir la date de début pour filtrer les deux premiers jours
    start_date = obs_temp.index.min() + pd.Timedelta(days=4)

    # Filtrer obs_temp pour supprimer les deux premiers jours
    obs_temp_filtered = obs_temp[obs_temp.index >= start_date]

    return sim_temp_filtered, obs_temp_filtered


def remove_first_days_sim(sim_temp, nb_delday):
    """
    Removes the first few days from the sim_temp DataFrame.

    Parameters:
    - sim_temp: DataFrame containing simulated temperature data with a 'Time' column in seconds.
    - nb_delday: Number of days to remove from the beginning of the simulation data.

    Returns:
    - sim_temp_filtered: Filtered sim_temp DataFrame.
    """
    # Filtrer sim_temp pour supprimer les deux premiers jours (86400 secondes par jour)
    sim_temp_filtered = sim_temp[sim_temp['Time'] >= 86400 * 4  ]

    return sim_temp_filtered

def reuse_end_in_inital():
    """
    Copies the 5th column of data from the 'S_pression_charge_temperature.dat' file 
    to the 'E_charge_initiale.dat' file.
    This function reads the source file, extracts the 5th column (index 4) of data, 
    and writes this column to the destination file.
    Raises:
        FileNotFoundError: If the source file does not exist.
    Comments:
        - Ensure that the source file exists before attempting to read it.
        - The 5th column is extracted using zero-based indexing (index 4).
        - The data is saved to the destination file in floating-point format.
    """
    
    
        # copy the colum 5 of the S_pression_charge_temperature.dat file in the E_charge_initiale.dat file
    source_file = '/home/ariviere/Programmes/ginette/application/2024_TD_ENS/SYNTHETIC_CASES/S_pression_charge_temperature.dat'
    destination_file = '/home/ariviere/Programmes/ginette/application/2024_TD_ENS/SYNTHETIC_CASES/E_charge_initiale.dat'

    # Read the source file
    source_data = np.loadtxt(source_file)
    if not os.path.exists(source_file):
        raise FileNotFoundError(f"{source_file} not found.")

    # Extract the 5th column (index 4)
    column_to_copy = source_data[:, 4]

    # Write the column to the destination file
    np.savetxt(destination_file, column_to_copy, fmt='%f')

