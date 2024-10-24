import pandas as pd
import numpy as np
import os
from pathlib import Path

def read_csv_with_multiple_separators(file_path, separators=[',', ';', '\t']):
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
            if df.shape[1] > 2:  # If more than one column, return the DataFrame
                return df
        except pd.errors.ParserError:
            continue
    raise ValueError(f"Cannot read the file {file_path} with the provided separators.")

import pandas as pd
import numpy as np

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
    
    # Extensive list of date formats
    formats = [
        # Common formats
        "%Y-%m-%d", "%d/%m/%Y", "%m/%d/%Y", "%Y/%m/%d",
        "%Y-%m-%d %H:%M:%S", "%d/%m/%Y %H:%M:%S", "%m/%d/%Y %H:%M:%S",
        "%Y/%m/%d %H:%M:%S", "%Y-%m-%d %H:%M", "%d/%m/%Y %H:%M", "%m/%d/%Y %H:%M", 
        "%Y/%m/%d %H:%M", "%Y-%m-%d %I:%M:%S %p", "%d/%m/%Y %I:%M:%S %p", "%m/%d/%Y %I:%M:%S %p",
        "%Y/%m/%d %I:%M:%S %p", "%Y-%m-%d %I:%M %p", "%d/%m/%Y %I:%M %p", "%m/%d/%Y %I:%M %p",
        "%Y/%m/%d %I:%M %p", "%d-%m-%Y", "%m-%d-%Y", "%Y.%m.%d", "%d.%m.%Y", "%m.%d.%Y", "%Y/%m/%d",
        
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
    
    # Try the generic method first (handles most common formats)
    try:
        df[date_column] = pd.to_datetime(df[date_column], errors='raise')
        return df
    except (ValueError, TypeError):
        pass  # If this fails, continue to try specific formats below

    # Iterate through the predefined formats
    for fmt in formats:
        try:
            df[date_column] = pd.to_datetime(df[date_column], format=fmt, errors='raise')
            return df  # Success, exit the function
        except (ValueError, TypeError):
            continue  # Try the next format
    
    # If no format works, raise a descriptive error
    raise ValueError(f"Failed to convert '{date_column}' to datetime with known formats.")

import os
from pathlib import Path
import pandas as pd

def process_obs_data(Obs_data, date_simul_bg, coef, ost, nb_day):
    all_data, temp_data = None, None
    
    # Reading the files and populating the dataframes
    for fichier in os.listdir(Obs_data):
        file_path = Path(Obs_data) / fichier 
        try:
            if 'deltaP' in fichier:
                all_data = read_csv_with_multiple_separators(file_path)
            elif 'Temp' in fichier:
                temp_data = read_csv_with_multiple_separators(file_path)
        except ValueError as e:
            print(f"Error processing file {fichier}: {e}")
    
    if all_data is None:
        raise ValueError("No 'deltaP' data could be read from the provided files.")
    
    if temp_data is not None:
        all_data = pd.merge(all_data, temp_data, on='dates', how='outer')
    
    # Drop unnecessary columns if they exist
    all_data = all_data.drop(columns=[col for col in ['#', 'ST2'] if col in all_data.columns], errors='ignore')
    
    # Convert 'dates' column to datetime
    all_data = convert_dates(all_data, 'dates')
    all_data = all_data[all_data['dates'] > date_simul_bg]

    # Ensure correct column names
    all_data.columns = ['dates', 'deltaP', 'TempMolo', 'Temp1', 'Temp2', 'Temp3', 'Temp4']

    # Reset index and process 'deltaP' values
    all_data.reset_index(drop=True, inplace=True)
    all_data['deltaP'] = all_data['deltaP'] * coef + ost

    # Define start and end dates
    date_begin = all_data['dates'].iloc[0]
    date_end = pd.to_datetime(date_begin) + pd.to_timedelta(nb_day, unit='d')
    print('Start date:', date_begin, 'End date:', date_end)

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
