# -*- coding: utf-8 -*-
"""
This module provides functions for reading, cleaning, merging, and interpolating observational data from CSV files,
with a focus on handling multiple date formats and aligning data to regular time intervals.
Functions:
----------
- read_csv_with_multiple_separators(file_path, separators=[',', ';', '\t']):
    Attempts to read a CSV file using a list of possible separators, returning a DataFrame if successful.
- convert_dates(df: pd.DataFrame, date_column: str) -> pd.DataFrame:
    Converts a specified column of string dates in a DataFrame to datetime objects, trying multiple common formats.
- process_obs_data(Obs_data, date_simul_bg, coef, ost, nb_day):
    Reads, merges, cleans, and interpolates observational data files (e.g., 'deltaP' and 'Temp'), aligns data to 15-minute intervals,
    applies scaling and offset to 'deltaP', and ensures regular time steps.
- process_obs_RIV2D(Station, Obs_data, date_simul_bg, nb_day, desc_station, pt100_coord):
    Processes and merges observational data for a given station, aligns data to 15-minute intervals, renames columns based on sensor
    mapping, and replaces temperature columns with indices from a provided mapping DataFrame.
------
- The module is designed to handle a variety of date formats and CSV file structures.
- Data is interpolated to fill missing values and ensure regular time intervals.
- Column renaming and mapping are performed to standardize sensor data across different files and stations.
"""
import pandas as pd
import numpy as np
import os
from pathlib import Path
import glob

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



def process_obs_data(Obs_data, date_simul_bg, coef, ost, nb_day):
    """
    Processes observational data by reading, merging, cleaning, and interpolating data files.
    Args:
        Obs_data (str): Path to the directory containing observational data files.
        date_simul_bg (datetime): The starting date for filtering the data.
        coef (float): Coefficient to scale the 'deltaP' values.
        ost (float): Offset to add to the scaled 'deltaP' values.
        nb_day (int): Number of days to calculate the end date from the start date.
    Returns:
        pd.DataFrame: A processed DataFrame with interpolated data aligned to 15-minute intervals.
    Raises:
        ValueError: If no 'deltaP' data could be read from the provided files.
    Notes:
        - The function reads files containing 'deltaP' and 'Temp' data, merges them, and processes the data.
        - Columns named '#' and 'ST2' are dropped if they exist.
        - The 'dates' column is converted to datetime and filtered based on `date_simul_bg`.
        - The 'deltaP' values are scaled and offset using `coef` and `ost`.
        - The data is reindexed to 15-minute intervals, and missing values are interpolated using time-based interpolation.
        - The function checks if all time differences between rows are 900 seconds (15 minutes) before and after interpolation.
    """

    all_data, temp_data = None, None
    
    # Reading the files and populating the dataframes
    for fichier in os.listdir(Obs_data):
        file_path = Path(Obs_data) / fichier 
        try:
            if 'deltaP' in fichier:
                all_data = read_csv_with_multiple_separators(file_path)
            elif 'Temp' in fichier:
                temp_data = read_csv_with_multiple_separators(file_path)
            elif 'Point' in fichier:
                all_data = read_csv_with_multiple_separators(file_path)
                print(f"Reading Point file: {fichier}")
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



def process_obs_RIV2D(Station, Obs_data, date_simul_bg, nb_day, desc_station, pt100_coord):
    """
    Processes observational data by reading, merging, cleaning, and interpolating data files.
    
    Args:
        Station (str): Name of the station.
        Obs_data (str): Path to the directory containing observational data files.
        date_simul_bg (datetime): Start date for filtering observations.
        nb_day (int): Number of days to calculate the end date from the start date.
        desc_station (dict): Dictionary mapping station codes to sensor names.
        pt100_coord (pd.DataFrame): DataFrame containing mapping of hobo and pt100 sensors to indices.

    Returns:
        pd.DataFrame: Cleaned and merged observational data.
    """
    # Construct the path to the station's data directory
    station_path = os.path.join(Obs_data, Station)
    all_data = None
    date_begin=date_simul_bg
    date_end = pd.to_datetime(date_begin) + pd.to_timedelta(nb_day, unit='d')
    print('Start date:', date_begin, 'End date:', date_end)
    # Iterate over all CSV files in the station's directory
    for csv_file in glob.glob(os.path.join(station_path, "*.csv")):
        file_name = os.path.basename(csv_file)
        # Extract the sensor type from the file name
        type_sensor = file_name.split("_")[0].replace(".csv", "")
        print("Processing sensor type:", type_sensor)

        # Read the CSV file using multiple possible separators
        df = read_csv_with_multiple_separators(csv_file, separators=[',', ';', '\t'])
        # Convert other columns than 'dates' to numeric
        for col in df.columns:
            if col != 'dates':
                df[col] = pd.to_numeric(df[col], errors='coerce')
        # Convert the 'dates' column to datetime format
        df = convert_dates(df, 'dates')
        # Ensure dates each have seconds set to zero 	2016-07-12 12:15:00 and not 	2016-07-12 12:15:02
        # Floor the dates to 15-minute intervals
        df['dates'] = df['dates'].dt.floor('15min')

        # Set 'dates' as the index for resampling
        df.set_index('dates', inplace=True)

        # Resample to 15-minute intervals and fill NaN values with hourly data if available
        df = df.resample('15T').mean()
        df.fillna(method='ffill', limit=3, inplace=True)  # Fill NaN using forward fill up to 1 hour (3 intervals)

        # Interpolate remaining missing values
        df.interpolate(method='time', inplace=True)

        # Reset the index to bring 'dates' back as a column
        df.reset_index(inplace=True)


        # Drop the conductivity column if it exists
        conductivity_column = next((col for col in df.columns if col.strip().lower() == 'conductivity [ms/cm]'), None)
        if conductivity_column:
            df.drop(columns=[conductivity_column], inplace=True)

        # Rename columns based on the sensor type
        if " Hydraulic Head [mNGF]" in df.columns:
            df.rename(columns={" Hydraulic Head [mNGF]": f"H_{type_sensor}"}, inplace=True)

        # Rename temperature depth columns
        for i in range(1, 5):
            col_name = f"temperature_depth_{i}_C"
            if col_name in df.columns:
                df.rename(columns={col_name: f"Temp_{type_sensor}_{i}"}, inplace=True)

        # Rename general temperature column
        if " Temperature [°C]" in df.columns:
            df.rename(columns={" Temperature [°C]": f"Temp_{type_sensor}"}, inplace=True)

        # Merge the current DataFrame with the accumulated data
        all_data = df if all_data is None else pd.merge(all_data, df, on='dates', how='outer')
    # find the sensor of RIV in desc_station
    if 'RIV' in desc_station:
        name_sensor_riv = desc_station['RIV']
        # create the column H_name_sensor_riv
        col_name_riv = f"H_{name_sensor_riv}"
    else:
        raise KeyError("'RIV' is not a key in desc_station.")


    hmax=all_data[col_name_riv].max()
    hmin=all_data[col_name_riv].mean()
    hqt95=all_data[col_name_riv].quantile(0.95)

    # Filter rows based on the start date
    all_data = all_data[all_data['dates'] >= date_simul_bg]
    all_data=all_data[all_data['dates'] <= date_end]

    # Replace sensor names with station codes using the desc_station dictionary
    for key, value in desc_station.items():
        all_data.columns = [
            col.replace(f"H_{value}", f"H_{key}").replace(f"Temp_{value}", f"T_{key}") if isinstance(col, str) else col
            for col in all_data.columns
        ]

    # Replace temperature columns in all_data with indices from pt100_coord
    pt100_coord["index"] = pt100_coord.index + 1  # Add an index column to pt100_coord
    for column in all_data.columns:
        if column.startswith("Temp_"):
            # Extract hobo and pt100 information from the column name
            parts = column.split("_")
            if len(parts) == 3 and parts[1].lower().startswith("hobo"):
                hobo = parts[1].lower()
                pt100 = int(parts[2])
                
                # Find the index in pt100_coord where hobo and pt100 match
                matching_row = pt100_coord[(pt100_coord['hobo'] == hobo) & (pt100_coord['pt100'] == pt100)]
                
                if not matching_row.empty:
                    # Extract the single index value
                    index_value = matching_row['index'].iloc[0]
                    # Rename the column with the corresponding index
                    all_data.rename(columns={column: f"Temp_{index_value}"}, inplace=True)

    # If 'H_RD' or 'H_RG' do not exist in all_data, copy the value from the other if available
    if 'H_RD' not in all_data.columns:
        all_data['H_RD'] = all_data['H_RG']
    if 'H_RG' not in all_data.columns:
        all_data['H_RG'] = all_data['H_RD']
    if 'T_RD' not in all_data.columns:
        all_data['T_RD'] = all_data['T_RG']
    if 'T_RG' not in all_data.columns:
        all_data['T_RG'] = all_data['T_RD']

    return all_data, hmax,hmin,hqt95