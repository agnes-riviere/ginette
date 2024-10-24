# Direct_model.py
import pandas as pd
import numpy as np
import subprocess

def setup_ginette(dt, state, nb_day, z_top, z_bottom, az, dz, date_simul_bg,dz_obs):
    print("la simulation commence à", date_simul_bg)
    
    # number of cell
    nb_cell=az/dz
    #-----------------------------------------------------------------
    ## write the setup of the modeled domain
    f_param_bck = open("E_parametre_backup.dat", "r")
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
    f_param_lec=open("E_cdt_aux_limites_bck.dat", "r")
    f_param_lec_new=open("E_cdt_aux_limites.dat", "w")
    lec_bc=f_param_lec.read()
    iclchgt = 1
    lec_bc = lec_bc.replace('[iclchgt]', '%1i' % iclchgt)
    lec_bc = lec_bc.replace('[itlecture]', '%08.0f' % dt)

    # Écrire les modifications dans le fichier
    f_param_lec_new.write(lec_bc)

    # Fermer les fichiers
    f_param_lec_new.close()
    f_param_lec.close()
    return z_obs
    
    #-----------------------------------------------------------------


def initial_conditions(all_data, z_top, z_bottom, dz, z_obs):
    # Initial conditions
    with open("E_temperature_initiale.dat", "w") as f_temp_IC:
        # Dynamically select temperature columns based on the names in all_data
        if 'T_top' in all_data.columns and 'T_bottom' in all_data.columns:
            temp_columns = ['T_top', 'T_bottom']
            z_temps = np.array([z_top + 0.005, z_bottom])
        else:
            temp_columns = ['TempMolo'] + [f'Temp{i+1}' for i in range(len(z_obs))]
            z_temps = np.array([z_top + 0.005] + z_obs)

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


    # To apply dp to the column, we need to convert it to a pressure value
    with open("E_charge_initiale.dat", "w") as f_chg_IC:
        initial_pres = all_data.iloc[0][['deltaP']] 
        initial_chg = pd.DataFrame({'z': [z_bottom], 'chg': [0]})
        new_row = pd.DataFrame({'z': [0], 'chg': initial_pres})
        initial_chg = pd.concat([initial_chg, new_row], ignore_index=True)
        initial_chg['z'] = initial_chg['z'].astype(float)
        initial_chg['chg'] = initial_chg['chg'].astype(float)
        z_values = np.arange(z_bottom + dz / 2, z_top + 0.001, dz)

        # Perform linear interpolation
        charge_init = np.interp(z_values, initial_chg['z'], initial_chg['chg'])  # Reversed to align with depths

        # Creating a DataFrame for the interpolated data
        interpolated_chg = pd.DataFrame({'z': z_values, 'chg': charge_init})
        interpolated_chg_sorted = interpolated_chg.sort_values(by='z', ascending=False)

        interpolated_chg_sorted['chg'].to_csv(f_chg_IC, index=False, sep='\n', header=False)



def boundary_conditions(all_data):
    # Boundary conditions
    all_data['bot'] = 0
    all_data['top'] = all_data['deltaP']
    
    # Save charge boundary conditions
    all_data[['top', 'bot']].to_csv('E_charge_t.dat', sep=' ', index=False, header=False)
    
    # Save temperature boundary conditions
    if 'T_top' in all_data.columns and 'T_bottom' in all_data.columns:
        all_data[['T_top', 'T_bottom']].to_csv('E_temp_t.dat', sep=' ', index=False, header=False)
    else:
        all_data[['TempMolo', 'Temp4']].to_csv('E_temp_t.dat', sep=' ', index=False, header=False)

# 

   
def generate_zone_parameters(z_bottom, dz, nb_zone, alt_thk, REF_k, REF_n, REF_l, REF_r, REF_k2=None, REF_n2=None, REF_l2=None, REF_r2=None):
    """
    Écrit le fichier de paramètres de zone en fonction du nombre de zones et des valeurs de paramètres.

    Parameters:
    - nb_zone: Nombre de zones (1 ou 2).
    - REF_k: Perméabilité intrinsèque pour la zone 1.
    - REF_n: Porosité pour la zone 1.
    - REF_l: Conductivité thermique pour la zone 1.
    - REF_r: Densité pour la zone 1.
    - REF_k2: Perméabilité intrinsèque pour la zone 2 (si nb_zone == 2).
    - REF_n2: Porosité pour la zone 2 (si nb_zone == 2).
    - REF_l2: Conductivité thermique pour la zone 2 (si nb_zone == 2).
    - REF_r2: Densité pour la zone 2 (si nb_zone == 2).
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
        # Lire le fichier de sauvegarde des paramètres de zone
        with open("E_zone_parameter_backup.dat", "r") as f_paramZ_bck:
            param_zone = f_paramZ_bck.read()
        # Remplacer les valeurs des paramètres
        param_zone = param_zone.replace('[k1]', '%8.2e' % REF_k)
        param_zone = param_zone.replace('[n1]', '%6.2f' % REF_n)
        param_zone = param_zone.replace('[l1]', '%6.2f' % REF_l)
        param_zone = param_zone.replace('[r1]', '%6.2f' % REF_r)
    elif nb_zone == 2:
        # Lire le fichier de sauvegarde des paramètres de zone pour 2 zones
        with open("E_zone_parameter_backup_2zones.dat", "r") as f_paramZ_bck:
            param_zone = f_paramZ_bck.read()
        # Remplacer les valeurs des paramètres pour les deux zones
        param_zone = param_zone.replace('[k1]', '%8.2e' % REF_k)
        param_zone = param_zone.replace('[n1]', '%6.2f' % REF_n)
        param_zone = param_zone.replace('[l1]', '%6.2f' % REF_l)
        param_zone = param_zone.replace('[r1]', '%6.2f' % REF_r)
        param_zone = param_zone.replace('[k2]', '%8.2e' % REF_k2)
        param_zone = param_zone.replace('[n2]', '%6.2f' % REF_n2)
        param_zone = param_zone.replace('[l2]', '%6.2f' % REF_l2)
        param_zone = param_zone.replace('[r2]', '%6.2f' % REF_r2)
    # Ouvrir le fichier de paramètres de zone
    f_paramZ_new = open("E_zone_parameter.dat", 'w')
    # Écrire les paramètres dans le nouveau fichier
    f_paramZ_new.write(param_zone)
    f_paramZ_new.close()

def run_direct_model(date_simul_bg,z_bottom, dz, nb_zone, alt_thk, REF_k, REF_n, REF_l, REF_r, REF_k2=None, REF_n2=None, REF_l2=None, REF_r2=None):
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



def remove_first_two_days(sim_temp, obs_temp):
    """
    Supprime les deux premiers jours des DataFrames sim_temp et obs_temp.

    Parameters:
    - sim_temp: DataFrame contenant les données de température simulées avec une colonne 'Time' en secondes.
    - obs_temp: DataFrame contenant les données de température observées avec un index de dates.

    Returns:
    - sim_temp_filtered: DataFrame sim_temp filtré.
    - obs_temp_filtered: DataFrame obs_temp filtré.
    """
    # Filtrer sim_temp pour supprimer les deux premiers jours (86400 secondes par jour)
    sim_temp_filtered = sim_temp[sim_temp['Time'] >= 86400 * 2]

    # Convertir l'index de obs_temp en format datetime si ce n'est pas déjà fait
    obs_temp.index = pd.to_datetime(obs_temp.index)

    # Définir la date de début pour filtrer les deux premiers jours
    start_date = obs_temp.index.min() + pd.Timedelta(days=2)

    # Filtrer obs_temp pour supprimer les deux premiers jours
    obs_temp_filtered = obs_temp[obs_temp.index >= start_date]

    return sim_temp_filtered, obs_temp_filtered

