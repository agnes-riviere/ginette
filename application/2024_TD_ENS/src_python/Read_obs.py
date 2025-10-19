import pandas as pd
import numpy as np
import os
from pathlib import Path

def read_csv_with_multiple_separators(file_path, separators=[';', ',', '\t']):
    """
    Try reading a CSV file with multiple separators until one works.
    
    Parameters:
    - file_path: Path to the CSV file.
    - separators: List of separators to try.
    
    Returns:
    - DataFrame: The DataFrame read from the CSV file.
    """
    for sep in separators:
        try:
            df = pd.read_csv(file_path, sep=sep)
            # Check if we have more than 1 column (not just everything in one column)
            if df.shape[1] > 1:  
                return df
        except (pd.errors.ParserError, UnicodeDecodeError):
            continue
    raise ValueError(f"Cannot read the file {file_path} with the provided separators.")

def convert_dates(df: pd.DataFrame, date_column: str) -> pd.DataFrame:
    """
    Convert dates from a list of strings in a specified column of the dataframe 
    by testing multiple common date formats. If none match, use the generic `pd.to_datetime` method.
    
    This function modifies the dataframe in place and raises an error if conversion fails.
    
    Parameters:
    -----------
    df : pd.DataFrame
        The input dataframe containing the date column.
    date_column : str
        The name of the column containing dates as strings.
    
    Returns:
    --------
    pd.DataFrame
        The dataframe with the date column converted to datetime objects.

    Raises:
    -------
    ValueError:
        If the date column cannot be converted with any known formats or the generic method.
    """
    
    # Validate input
    if df is None or date_column not in df.columns:
        raise ValueError(f"DataFrame is None or column '{date_column}' does not exist.")
    
    # Debug: Print detailed information about the dates
    print(f"Sample dates to convert: {df[date_column].head().tolist()}")
    print(f"Date column dtype: {df[date_column].dtype}")
    
    # Check for any null values
    null_count = df[date_column].isnull().sum()
    if null_count > 0:
        print(f"Warning: Found {null_count} null values in date column")
    
    # Try the mixed format approach first (handles inconsistent formats within same column)
    try:
        df[date_column] = pd.to_datetime(df[date_column], format='mixed', dayfirst=True)
        print("✓ Successfully converted dates using format='mixed' with dayfirst=True")
        return df
    except (ValueError, TypeError) as e:
        print(f"Mixed format conversion failed: {e}")
    
    # Try the generic method with dayfirst=True
    try:
        df[date_column] = pd.to_datetime(df[date_column], errors='raise', dayfirst=True)
        print("✓ Successfully converted dates using generic pd.to_datetime with dayfirst=True")
        return df
    except (ValueError, TypeError) as e:
        print(f"Generic conversion failed: {e}")

    # If mixed and generic methods fail, try specific formats
    priority_formats = ["%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%d/%m/%Y"]
    
    for fmt in priority_formats:
        try:
            print(f"Trying format: {fmt}")
            df[date_column] = pd.to_datetime(df[date_column], format=fmt, errors='raise')
            print(f"✓ Successfully converted dates using format: {fmt}")
            return df
        except (ValueError, TypeError) as e:
            print(f"Format {fmt} failed: {e}")
            continue

    # If priority formats fail, try all other formats
    other_formats = [
        # Common formats
        "%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d",
        "%Y-%m-%d %H:%M:%S", "%m/%d/%Y %H:%M:%S",
        "%Y/%m/%d %H:%M:%S", "%Y-%m-%d %H:%M", "%m/%d/%Y %H:%M", 
        "%Y/%m/%d %H:%M", "%Y-%m-%d %I:%M:%S %p", "%d/%m/%Y %I:%M:%S %p", "%m/%d/%Y %I:%M:%S %p",
        "%Y/%m/%d %I:%M:%S %p", "%Y-%m-%d %I:%M %p", "%d/%m/%Y %I:%M %p", "%m/%d/%Y %I:%M %p",
        "%Y/%m/%d %I:%M %p",
        
        # Year-month formats
        "%Y-%m", "%Y/%m", "%m/%Y", "%m-%Y", "%Y.%m",
        
        # Short year variations
        "%d/%m/%y", "%m/%d/%y", "%y/%m/%d", "%y-%m-%d", "%y.%m.%d", "%d-%m-%y", "%m-%d-%y",
        "%d/%b/%y", "%d-%b-%y", "%d %b %y", "%d-%B-%y", "%d/%B/%y",
        
        # Full month names
        "%d %B %Y", "%d %b %Y", "%B %d, %Y", "%b %d, %Y", "%B %d %Y", "%b %d %Y",
        
        # Time only formats
        "%H:%M:%S", "%I:%M:%S %p", "%H:%M", "%I:%M %p",
        
        # Mixed separators
        "%Y.%m.%d %H:%M:%S", "%d.%m.%Y %H:%M:%S", "%m.%d.%Y %H:%M:%S",
        "%Y.%m.%d %I:%M:%S %p", "%d.%m.%Y %I:%M:%S %p", "%m.%d.%Y %I:%M:%S %p",
        "%Y.%m.%d %H:%M", "%d.%m.%Y %H:%M", "%m.%d.%Y %H:%M", "%Y/%m/%d %I:%M %p",
        
        # Variations of month-day-year
        "%m-%d-%Y %H:%M:%S", "%m-%d-%Y %I:%M:%S %p", "%m-%d-%Y %H:%M",
        "%m-%d-%Y %I:%M %p", "%m/%d/%Y %I:%M %p", "%m.%d.%Y %I:%M %p",
        
        # Day first vs Month first ambiguities
        "%d-%b-%Y", "%b-%d-%Y", "%d-%B-%Y", "%B-%d-%Y",
        
        # ISO formats
        "%Y%m%dT%H%M%S", "%Y%m%dT%H%M", "%Y%m%d %H%M%S", "%Y%m%d %H%M", "%Y%m%d",
        
        # Others (edge cases)
        "%A, %B %d, %Y", "%a, %b %d, %Y", "%A, %B %d %Y", "%a, %b %d %Y",
        "%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S",
        "%Y/%m/%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S.%f"
    ]
    
    for fmt in other_formats:
        try:
            df[date_column] = pd.to_datetime(df[date_column], format=fmt, errors='raise')
            print(f"✓ Successfully converted dates using format: {fmt}")
            return df
        except (ValueError, TypeError):
            continue
    
    # If no format works, raise a descriptive error with more information
    sample_dates = df[date_column].head().tolist()
    raise ValueError(f"Failed to convert '{date_column}' to datetime with known formats. Sample dates: {sample_dates}")

def process_obs_data(Obs_data, date_simul_bg, coef, offset, nb_day):  # Fixed parameter name: ost -> offset
    """
    Traite les données d'observation spécifiques au Point1Touques
    
    Parameters:
    -----------
    Obs_data : str
        Répertoire contenant les fichiers CSV
    date_simul_bg : pd.Timestamp  # Fixed parameter name
        Date de début de la simulation
    coef : float
        Coefficient multiplicateur pour la pression
    offset : float  # Fixed parameter name
        Décalage pour la pression
    nb_day : int  # Fixed parameter name
        Durée de simulation en jours
    
    Returns:
    --------
    pd.DataFrame
        DataFrame avec index temporel et colonnes Temp1, Temp2, Temp3, Pressure
    """
    all_data, temp_data = None, None
    
    # Reading the files and populating the dataframes
    for fichier in os.listdir(Obs_data):
        file_path = Path(Obs_data) / fichier 
        try:
            if 'deltaP' in fichier:
                all_data = read_csv_with_multiple_separators(file_path)
                print(f"Processing file: {file_path}")
                print(all_data.head())
                # Clean up the deltaP data - remove the first column if it contains row numbers
                if all_data.columns[0] in ['#', 'Unnamed: 0'] or all_data.iloc[0, 0] == '1':
                    all_data = all_data.drop(all_data.columns[0], axis=1)
                # Rename columns to standard names
                expected_cols = ['dates', 'deltaP', 'T']
                if len(all_data.columns) >= len(expected_cols):
                    all_data.columns = expected_cols[:len(all_data.columns)]
                
            elif 'Temp' in fichier:
                temp_data = read_csv_with_multiple_separators(file_path)
                print(f"Processing file: {file_path}")
                print(temp_data.head())
                # Clean up the temperature data - remove the first column if it contains row numbers
                if temp_data.columns[0] in ['ST2', 'Unnamed: 0'] or temp_data.iloc[0, 0] == '1':
                    temp_data = temp_data.drop(temp_data.columns[0], axis=1)
                # Rename columns to standard names
                expected_cols = ['dates', 'T1', 'T2', 'T3', 'T4']
                if len(temp_data.columns) >= len(expected_cols):
                    temp_data.columns = expected_cols[:len(temp_data.columns)]
                    
        except ValueError as e:
            print(f"Error processing file {fichier}: {e}")
    all_data=convert_dates(all_data, 'dates')
    temp_data=convert_dates(temp_data, 'dates')
    if all_data is None:
        raise ValueError("No 'deltaP' data could be read from the provided files.")
    
    if temp_data is not None:
        all_data = pd.merge(all_data, temp_data, on='dates', how='outer')
    
    # Print column names for debugging
    print("Noms des colonnes après lecture des fichiers:", list(all_data.columns))
    
    # Debug: Show sample of merged data before date conversion
    print("Sample of merged data before date conversion:")
    print(all_data[['dates']].head())
    
    # Convert 'dates' column to datetime
    all_data = convert_dates(all_data, 'dates')
    all_data = all_data[all_data['dates'] > date_simul_bg]

    # Ensure correct column names based on what we actually have
    if 'T1' in all_data.columns:
        # We have temperature data merged
        all_data = all_data.rename(columns={
            'T': 'TempMolo',
            'T1': 'Temp1', 
            'T2': 'Temp2', 
            'T3': 'Temp3', 
            'T4': 'Temp4'
        })
    else:
        # Only pressure data, rename T column
        all_data = all_data.rename(columns={'T': 'TempMolo'})
        # Add missing temperature columns with NaN
        for col in ['Temp1', 'Temp2', 'Temp3', 'Temp4']:
            if col not in all_data.columns:
                all_data[col] = np.nan

    # Reset index and process 'deltaP' values
    all_data.reset_index(drop=True, inplace=True)
    all_data['deltaP'] = all_data['deltaP'] * coef + offset  # Fixed variable name: ost -> offset

    # Define start and end dates
    date_begin = all_data['dates'].iloc[0]
    date_end = pd.to_datetime(date_begin) + pd.to_timedelta(nb_day, unit='d')  # Fixed variable name

    # Time difference calculation
    time_diff = all_data['dates'].diff().dropna()
    indices_not_equal_900_sec = time_diff[time_diff != pd.Timedelta(seconds=900)].index

    # Check if all differences are 900 seconds (15 minutes)
    is_equal_900_sec = all(time_diff == pd.Timedelta(seconds=900))

    # Floor the dates to 15-minute intervals
    all_data['dates'] = all_data['dates'].dt.floor('15min')

    # Set 'timestamps' column for reindexing and interpolation
    all_data['timestamps'] = all_data['dates']

    # Create a new index with 15-minute intervals
    new_index = pd.date_range(start=all_data['timestamps'].min(), end=all_data['timestamps'].max(), freq='15min')

    # Reindex DataFrame to align with new 15-minute intervals
    all_data = all_data.set_index('timestamps').reindex(new_index)

    # Ensure 'timestamps' as index
    all_data.index.name = 'timestamps'

    # Columns for interpolation
    columns_to_interpolate = ['deltaP', 'TempMolo', 'Temp1', 'Temp2', 'Temp3', 'Temp4']

    # Apply time-based interpolation
    all_data[columns_to_interpolate] = all_data[columns_to_interpolate].interpolate(method='time')

    # Time difference check after interpolation
    time_diff_after_interp = all_data.index.to_series().diff().dropna()
    indices_not_equal_900_sec_after_interp = time_diff_after_interp[time_diff_after_interp != pd.Timedelta(seconds=900)].index

    # Results display
    if len(indices_not_equal_900_sec_after_interp) == 0:
        print("Toutes les différences entre les lignes sont de 900 secondes après interpolation.")
    else:
        print("Les différences ne sont pas toutes de 900 secondes après interpolation. Voici les indices concernés :")
        print(indices_not_equal_900_sec_after_interp)
        print(all_data.loc[indices_not_equal_900_sec_after_interp])
    
    return all_data
