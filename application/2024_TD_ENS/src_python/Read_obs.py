import os
import pandas as pd
import numpy as np

def rename_columns(df):
    """
    Rename columns based on specific conditions.
    
    Parameters:
    - df: pandas DataFrame with columns to be renamed.
    
    Returns:
    - df: pandas DataFrame with renamed columns.
    """
    def column_mapper(col_name):
        if 'temperature1' in col_name.lower():
            return 'Temp1'
        elif 'temperature2' in col_name.lower():
            return 'Temp2'
        elif 'stream' in col_name.lower() or 'molo' in col_name.lower():
            return 'TempMolo'
        else:
            return col_name

    df = df.rename(columns=column_mapper)
    return df


def convert_dates(df: pd.DataFrame, date_column: str):
    """
    Convert dates from a list of strings by testing several different input formats.
    Try all date formats already encountered in data points.
    If none of them is OK, try the generic way (None).
    If the generic way doesn't work, this method fails
    (in that case, you should add the new format to the list).
    
    This function works directly on the given Pandas dataframe (in place).
    This function assumes that the specified column of the given Pandas dataframe
    contains the dates as character string type.
    
    For datetime conversion performance, see:
    See https://stackoverflow.com/questions/40881876/python-pandas-convert-datetime-to-timestamp-effectively-through-dt-accessor
    """
    formats = [
        "%m/%d/%y %H:%M:%S", "%m/%d/%y %I:%M:%S %p",
        "%d/%m/%y %H:%M",    "%d/%m/%y %I:%M %p",
        "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %I:%M:%S %p", 
        "%d/%m/%Y %H:%M",    "%d/%m/%Y %I:%M %p",
        "%y/%m/%d %H:%M:%S", "%y/%m/%d %I:%M:%S %p", 
        "%y/%m/%d %H:%M",    "%y/%m/%d %I:%M %p",
        "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %I:%M:%S %p", 
        "%Y/%m/%d %H:%M",    "%Y/%m/%d %I:%M %p",
        "%d/%m/%Y %H:%M",    "%d/%m/%Y %H:%M:%S", "%Y-%m-%d %H:%M:%S",
        None
    ]
    
    times = df[date_column]
    for f in formats:
        try:
            # Convert strings to datetime objects
            new_times = pd.to_datetime(times, format=f)
            # Convert datetime series to numpy array of integers (timestamps)
            new_ts = new_times.values.astype(np.int64)
            # If times are not ordered, this is not the appropriate format
            test = np.sort(new_ts) - new_ts
            if np.sum(abs(test)) != 0:
                raise ValueError()
            # Else, the conversion is a success
            df[date_column] = new_times
            return
        
        except (ValueError, TypeError):
            continue
    
    # None of the known formats are valid
    raise ValueError("Cannot convert dates: No known formats match your data!")
    return df


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
    # Convert dates in the 'dates' column
    convert_dates(all_data, 'dates')
    all_data = all_data[all_data['dates'] > date_simul_bg]
    convert_dates(temp_data,'dates')
    temp_data = temp_data[temp_data['dates'] > date_simul_bg]

    all_data = all_data.merge(temp_data, on=['dates'])
    # drop 2 columns # and ST2
    all_data = all_data.drop(columns=['ST2', '#'], inplace=False)

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

