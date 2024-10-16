import os
import pandas as pd

def process_obs_data(Obs_data, date_simul_bg, coef, ost, nb_day):
    # Obtenez la liste des fichiers dans le répertoire
    fichiers = os.listdir(Obs_data)
    print(fichiers)
    info = None
    all_data = None

    for fichier in fichiers:
        if fichier.endswith('.csv'):
            if 'info' in fichier:
                info = pd.read_csv(os.path.join(Obs_data, fichier))
            elif 'deltaP' in fichier:
                all_data = pd.read_csv(os.path.join(Obs_data, fichier), sep=';')
            elif 'Temp' in fichier:
                temp_data = pd.read_csv(os.path.join(Obs_data, fichier), sep=';')    

    all_data['dates'] = pd.to_datetime(all_data['dates'], format='%d/%m/%Y %H:%M')
    all_data = all_data[all_data['dates'] > date_simul_bg]

    temp_data['dates'] = pd.to_datetime(temp_data['dates'], format='%d/%m/%Y %H:%M:%S')
    temp_data = temp_data[temp_data['dates'] > date_simul_bg]

    all_data = all_data.merge(temp_data, on=['dates'])
    # drop 2 columns # and ST2
    all_data = all_data.drop(columns=['ST2', '#'], inplace=False)
    # rename columns
    all_data.columns = ['dates', 'deltaP', 'TempMolo', 'Temp1', 'Temp2', 'Temp3', 'Temp4']

    # Réindexer le DataFrame
    all_data.reset_index(drop=True, inplace=True)
    all_data['deltaP'] = all_data['deltaP'] * coef + ost

    date_begin = all_data['dates'].iloc[0]
    date_end = pd.to_datetime(date_begin) + pd.to_timedelta(nb_day, unit='d')
    print(date_begin, date_end)
    # pas de temps
    time_diff = all_data['dates'].diff().dropna()
    indices_not_equal_900_sec = time_diff[time_diff != pd.Timedelta(seconds=900)].index

    # pas de temps=900s
    is_equal_900_sec = all(time_diff == pd.Timedelta(seconds=900))
    
    # Trouver les horodatages qui ne commencent pas à des intervalles de 15 minutes et ajuster uniquement ceux-là
#mask = all_data['dates'].dt.minute % 15 != 0
#all_data.loc[mask, 'timestamps'] = all_data.loc[mask, 'dates'] - pd.to_timedelta(all_data.loc[mask, 'dates'].dt.minute % 15, unit='m')
# Réinitialiser les secondes à zéro pour tous les horodatages
    all_data['dates'] = all_data['dates'].dt.floor('min')
    all_data['dates'] = all_data['dates'].dt.floor('15T')
    all_data['timestamps']=all_data['dates']

# Créez un nouvel index avec des dates toutes les 15 minutes
    new_index = pd.date_range(start=all_data['timestamps'].min(), end=all_data['timestamps'].max(), freq='15T')

# Réindexez le DataFrame pour inclure ces nouvelles dates
    all_data = all_data.set_index('timestamps').reindex(new_index)
# Sélection des colonnes pour l'interpolation
    columns_to_interpolate = ['deltaP', 'TempMolo', 'Temp1', 'Temp2', 'Temp3', 'Temp4']

# Appliquez l'interpolation uniquement sur la colonne 'valeur'
    all_data[columns_to_interpolate] = all_data[columns_to_interpolate].interpolate(method='time')
    
# Vérification des écarts temporels après l'interpolation
    time_diff = all_data['dates'].diff().dropna()
    indices_not_equal_900_sec = time_diff[time_diff != pd.Timedelta(seconds=900)].index


    if len(indices_not_equal_900_sec) == 0:
        print("Toutes les différences entre les lignes sont de 900 secondes.")
    else:
        print("Les différences ne sont pas toutes de 900 secondes. Voici les indices concernés :")
        print(indices_not_equal_900_sec)
        # Afficher les lignes correspondantes dans le dataframe
        print(all_data.loc[indices_not_equal_900_sec])
    
    return all_data

