import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm

def plot_obs(all_data):
    """
    Plots observation data for 'dp' and multiple temperature readings.
    This function creates a figure with two vertically stacked subplots:
    one for 'dp' (delta pressure) and one for various temperature readings.
    The 'dp' subplot displays a single line plot, while the temperature subplot
    displays multiple line plots for different temperature columns.
    Parameters:
    all_data (DataFrame): A pandas DataFrame containing the data to be plotted.
                        It must include the following columns:
                        - 'dates': The dates for the x-axis.
                        - 'deltaP': The delta pressure values for the 'dp' plot.
                        - 'TempMolo', 'Temp1', 'Temp2', 'Temp3', 'Temp4': The temperature values for the temperature plot.
    Returns:
    None: This function does not return any value. It displays the plots.
    """
    
   
    # Creating subplots for 'dp' and 'temp' graphs vertically
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Plotting 'dp'
    ax1.plot(all_data['dates'], all_data['deltaP'], 'b-', label='dp')
    ax1.set_ylabel('dp')
    ax1.set_title('dp')
    ax1.legend()

    # Plotting 'Temps'
    temp_columns = ['TempMolo', 'Temp1', 'Temp2', 'Temp3', 'Temp4']
    colors = ['r', 'g', 'b', 'c', 'm']
    for i, col in enumerate(temp_columns):
        ax2.plot(all_data['dates'], all_data[col], color=colors[i], label=col)

    ax2.set_xlabel('Date')
    ax2.set_ylabel('Temperature')
    ax2.set_title('Temperature')
    ax2.legend()

    # Adjust layout to prevent overlap
    plt.tight_layout()

    plt.show()


def plot_obs_zoom(all_data, start_date, end_date):
    """
    Plots zoomed-in observations of 'dp' and 'Temperature' within a specified date range.
    This function filters the input data based on the provided start and end dates,
    and then creates two subplots: one for 'dp' and another for various temperature readings.
    Parameters:
    all_data (DataFrame): The input data containing 'dates', 'deltaP', and temperature columns.
    start_date (str or datetime): The start date for filtering the data.
    end_date (str or datetime): The end date for filtering the data.
    Returns:
    None: This function does not return any value. It displays the plots.
    """
    
    # Filtrer les données selon les dates spécifiées
    filtered_data = all_data[(all_data['dates'] > start_date) & (all_data['dates'] < end_date)]
    
    # Creating subplots for 'dp' and 'temp' graphs vertically
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Plotting 'dp'
    ax1.plot(filtered_data['dates'], filtered_data['deltaP'], 'b-', label='dp')
    ax1.set_ylabel('dp')
    ax1.set_title('dp')
    ax1.legend()

    # Plotting 'Temps'
    temp_columns = ['TempMolo', 'Temp1', 'Temp2', 'Temp3', 'Temp4']
    colors = ['r', 'g', 'b', 'c', 'm']
    for i, col in enumerate(temp_columns):
        ax2.plot(filtered_data['dates'], filtered_data[col], color=colors[i], label=col)

    ax2.set_xlabel('Date')
    ax2.set_ylabel('Temperature')
    ax2.set_title('Temperature')
    ax2.legend()

    # Adjust layout to prevent overlap
    plt.tight_layout()

    plt.show()

def plot_domain(nb_zone, alt_thk, z_top, z_bottom):
    """
    Plot cells based on coordinates x and z, and color the cells based on the condition alt_thk > z.

    Parameters:
    - nb_zone (int): Number of zones. If nb_zone=1, the porous medium is homogeneous.
    - alt_thk (float): Altitude threshold for coloring the cells.
    - z_top (float): Top z-coordinate for the plot.
    - z_bottom (float): Bottom z-coordinate for the plot.
    """
    # Read coordinates with z_top and z_bottom
    coord = pd.read_csv('E_coordonnee.dat', sep='\s+', header=None, names=['x', 'z'])
    zone_parameters = pd.read_csv('E_zone_parameter.dat', sep='\s+', header=None, names=['zone', 'k', 'n', 'l', 'r'])

    # Extract coordinates
    x = coord['x'].values
    z = coord['z'].values

    # Calculate cell height
    z_unique = np.sort(np.unique(z))
    cell_height = np.diff(z_unique).mean() if len(z_unique) > 1 else 1.0  # Use 1.0 by default if only one value

    # Define cell width (you can adjust this value)
    x_unique = np.sort(np.unique(x))
    cell_width = np.diff(x_unique).mean() if len(x_unique) > 1 else 1.0  # Use 1.0 by default if only one value

    # Create a grid of cells centered on the points (x, z)
    fig, ax = plt.subplots(figsize=(10, 6))
    # Plot cells with colors based on the condition alt_thk > z or alt_thk >= z


    # Plot cells with colors based on the condition alt_thk > z
    for xi, zi in zip(x, z):
        if nb_zone == 1:
            color = 'blue'  # Color for homogeneous porous medium
        else:
            if zi >= alt_thk:
                color = 'red'  # Color if alt_thk > z
            else:
                color = 'green'  # Color otherwise
        rect = plt.Rectangle((xi - cell_width / 2, zi - cell_height / 2), cell_width, cell_height, edgecolor='black', facecolor=color)
        ax.add_patch(rect)

    # Add labels to the points
    for i, (xi, zi) in enumerate(zip(x, z)):
        ax.text(xi, zi, str(i), ha='center', va='center', color='white')

    # Configure axes
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    ax.set_title('Domain')
    ax.set_aspect(10)  # Adjust aspect ratio so that z scale is 10 times x scale
    ax.set_ylim(z_bottom, z_top)  # Adjust y-axis to go from z_bottom to z_top
    plt.grid(True)
    plt.show()


def plot_compare_temperatures_obs_sim(sim_temp, obs_temp,fontsize=15):
    """
    Plots and compares observed and simulated temperatures for each sensor.

    Parameters:
    - sim_temp: DataFrame containing the simulated temperatures. Must contain a 'dates' column and temperature columns named 'Temp1', 'Temp2', 'Temp3', etc.
    - obs_temp: DataFrame containing the observed temperatures. Must contain a 'dates' column and temperature columns named 'Temp1', 'Temp2', 'Temp3', etc.
    - fontsize: Font size for the labels (default is 15).
    """
    zoomSize = 2
    titleSize = fontsize + zoomSize
    fig, axes = plt.subplots(1, 3, figsize=(20, 5), sharey=True)

    # Assurez-vous que les dates sont alignées
    sim_temp['dates'] = pd.to_datetime(sim_temp['dates'])
    obs_temp['dates'] = pd.to_datetime(obs_temp['dates'])
    sim_temp.set_index('dates', inplace=True)
    obs_temp.set_index('dates', inplace=True)

    # Fusionner les DataFrames pour s'assurer qu'ils ont la même taille
    merged_df = pd.merge(sim_temp, obs_temp, left_index=True, right_index=True, suffixes=('_sim', '_obs'))

    # Colonnes de température (en supposant qu'elles soient nommées 'Temp1', 'Temp2', 'Temp3')
    temp_columns = ['Temp1', 'Temp2', 'Temp3']

    axes[0].set_ylabel("Temperature in C", fontsize=fontsize)  # Label y-axis

    for i, col in enumerate(temp_columns):
        axes[i].set_xlabel("Date", fontsize=fontsize)
        axes[i].plot(merged_df.index, merged_df[f'{col}_sim'], label="Simulated")
        axes[i].plot(merged_df.index, merged_df[f'{col}_obs'], label="Observed")
        axes[i].legend()
        axes[i].set_title(f"Sensor {i+1}")

    plt.subplots_adjust(wspace=0.05)
    plt.tight_layout()
    plt.show()

def plot_temperatures_sim(sim_temp, fontsize=15):
    """
    Plots the simulated temperatures for each sensor.
    Parameters:
    - sim_temp: DataFrame containing the simulated temperatures with a 'Date' column.
    - fontsize: Font size for the labels (default is 15).
    """
    zoomSize = 2
    titleSize = fontsize + zoomSize
    fig, axes = plt.subplots(1, 3, figsize=(20, 5), sharey=True)

    # Assurez-vous que les dates sont alignées
    sim_temp['dates'] = pd.to_datetime(sim_temp['dates'])
    sim_temp.set_index('dates', inplace=True)

    # Colonnes de température (en supposant qu'elles soient nommées 'Temp1', 'Temp2', 'Temp3')
    temp_columns = ['Temp1', 'Temp2', 'Temp3']

    axes[0].set_ylabel("Temperature in C", fontsize=fontsize)  # Label y-axis

    for i, col in enumerate(temp_columns):
        axes[i].set_xlabel("Date", fontsize=fontsize)
        axes[i].plot(sim_temp.index, sim_temp[f'{col}'], label="Simulated")
        axes[i].legend()
        axes[i].set_title(f"Sensor {i+1}")

    plt.subplots_adjust(wspace=0.05)
    plt.tight_layout()
    plt.show()

def plot_temperatures_profiles(fontsize=15, interval=43200):
    """
    This function plots the evolution of the temperature profile over time at specified intervals.
    Plots the evolution of the temperature profile at specified intervals.

    Parameters:
    - interval (int): Time interval in seconds for displaying the temperature profiles (default is 43200 seconds, or 12 hours).
    - interval: Time interval in seconds for displaying the temperature profiles (default is 43200 seconds, or 12 hours).
    """
    sim_profile = 'Sim_temperature_profil_t.dat'
    # Lire les profils de température
    sim_temp = pd.read_csv(sim_profile, sep='\s+', header=None, names=['Time', 'z', 'Temp'])

    # Filtrer les données pour n'imprimer que toutes les 'interval' secondes
    filtered_sim_temp = sim_temp[sim_temp['Time'] % interval == 0]

    # Plot de l'évolution du profil de température à chaque intervalle spécifié
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    for time in filtered_sim_temp['Time'].unique():
        temp_profile = filtered_sim_temp[filtered_sim_temp['Time'] == time]
        if interval >= 86400:  # Si l'intervalle dépasse un jour
            ax.plot(temp_profile['Temp'], temp_profile['z'], label=f'Time = {time/86400:.1f} days')
        else:
            ax.plot(temp_profile['Temp'], temp_profile['z'], label=f'Time = {time/3600:.1f} h')

    ax.set_ylabel("Depth (m)", fontsize=fontsize)
    ax.set_xlabel("T (°C)", fontsize=fontsize)
    ax.grid()
    if interval >= 86400:
        ax.set_title(f"Evolution du profil de température tous les {interval/86400:.1f} jours", fontsize=fontsize, pad=20)
    else:
        ax.set_title(f"Evolution du profil de température toutes les {interval/3600:.1f} heures", fontsize=fontsize, pad=20)
    # Placer la légende à l'extérieur de la figure
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=fontsize)
    plt.tight_layout()
    plt.show()


def plot_water_profiles_interpol(fontsize=15, date_simul_bg=None):
    """
    The above code is a Python docstring that provides a description of a function that plots the water
    velocity profile as a function of depth (z) and dates. It specifies the parameters `fontsize` and
    `date_simul_bg`, with default values for `fontsize` set to 15. The `date_simul_bg` parameter
    represents the start date of the simulation.
    Plots the water velocity profile as a function of depth (z) and dates.

    Parameters:
    - fontsize: Font size for the labels (default is 15).
    - date_simul_bg (str or datetime, optional): Start date of the simulation in a format recognized by pandas.to_datetime. If not provided, defaults to the current date and time.
    """
    # Set default date_simul_bg if not provided
    if date_simul_bg is None:
        date_simul_bg = pd.to_datetime('now')
    else:
        date_simul_bg = pd.to_datetime(date_simul_bg)

    # Lire les profils de vitesse
    sim_water_discharge = pd.read_csv('Sim_velocity_profil_t.dat', sep='\s+', header=None, names=['Time', 'z', 'q'])
    # delete first day of simulation
    sim_water_discharge = sim_water_discharge[sim_water_discharge['Time'] > 86400*2]
    sim_water_discharge['dates'] = pd.to_datetime(sim_water_discharge['Time'], unit='s', origin=date_simul_bg)

    # Supprimer les doublons
    sim_water_discharge = sim_water_discharge.drop_duplicates(subset=['q', 'z', 'Time'])

    # Créer une grille pour l'interpolation
    z_new = np.linspace(sim_water_discharge['z'].min(), sim_water_discharge['z'].max(), num=500)
    t_new = np.linspace(sim_water_discharge['Time'].min(), sim_water_discharge['Time'].max(), num=500)
    z_grid, t_grid = np.meshgrid(z_new, t_new)

    # Interpolation des données
    points = np.array([sim_water_discharge['z'], sim_water_discharge['Time']]).T
    values = sim_water_discharge['q']
    q_interp = griddata(points, values, (z_grid, t_grid), method='cubic')
    # t_grid in datetime format
    t_new = pd.to_datetime(t_new, unit='s', origin=date_simul_bg)
    # Tracer les profils interpolés
    fig, ax = plt.subplots(1, 1, figsize=(15, 8))
    c = ax.contourf(t_new, z_new, q_interp.T, cmap='Spectral_r', levels=100)
    fig.colorbar(c, ax=ax, label='Water velocity (m/s)')

    ax.set_ylabel("Depth (m)", fontsize=fontsize)
    ax.set_xlabel("Date", fontsize=fontsize)
    ax.grid()
    ax.set_title("Profil de vitesse d'eau", fontsize=fontsize, pad=20)
    
    # Formater les dates sur l'axe des x
    ax.xaxis_date()
    fig.autofmt_xdate()

    plt.tight_layout()
    plt.show()

def plot_temperature_profiles_interpol(fontsize=15, date_simul_bg=None):
    """
    This function plots the temperature profile over time and depth using interpolated data.
    Plots the temperature profile as a function of depth (z) and dates.

    Parameters:
    - fontsize: Font size for the labels (default is 15).
    - date_simul_bg: Start date of the simulation.
    """
    # Lire les profils de température
    sim_temp = pd.read_csv('Sim_temperature_profil_t.dat', sep='\s+', header=None, names=['Time', 'z', 'Temp'])
    # Supprimer les deux premiers jours de simulation
    sim_temp = sim_temp[sim_temp['Time'] > 86400 * 2]
    sim_temp['dates'] = pd.to_datetime(sim_temp['Time'], unit='s', origin=date_simul_bg)

    # Supprimer les doublons
    sim_temp = sim_temp.drop_duplicates(subset=['Temp', 'z', 'Time'])

    # Créer une grille pour l'interpolation
    z_new = np.linspace(sim_temp['z'].min(), sim_temp['z'].max(), num=500)
    t_new = np.linspace(sim_temp['Time'].min(), sim_temp['Time'].max(), num=500)
    z_grid, t_grid = np.meshgrid(z_new, t_new)

    # Interpolation des données
    points = np.array([sim_temp['z'], sim_temp['Time']]).T
    values = sim_temp['Temp']
    temp_interp = griddata(points, values, (z_grid, t_grid), method='cubic')
    # Convertir t_new en format datetime
    t_new = pd.to_datetime(t_new, unit='s', origin=date_simul_bg)
    
    # Tracer les profils interpolés
    fig, ax = plt.subplots(1, 1, figsize=(15, 8))
    c = ax.contourf(t_new, z_new, temp_interp.T, cmap='Spectral_r', levels=100)
    fig.colorbar(c, ax=ax, label='Temperature (°C)')

    ax.set_ylabel("Depth (m)", fontsize=fontsize)
    ax.set_xlabel("Date", fontsize=fontsize)
    ax.grid()
    ax.set_title("Profil de température", fontsize=fontsize, pad=20)
    
    # Formater les dates sur l'axe des x
    ax.xaxis_date()
    fig.autofmt_xdate()

    plt.tight_layout()
    plt.show()


def plot_heat_flux_profiles_interpolated(fontsize=15, date_simul_bg=None):
    """
    Plots the heat flux profiles (conductive, advective, and total) as a function of depth (z) and dates.
    
    Returns:
    None: This function does not return any value. It displays the plots.

    Parameters:
    - fontsize: Font size for the labels (default is 15).
    - date_simul_bg: Start date of the simulation.
    """
    # Lire les profils de flux de chaleur
    sim_flux = pd.read_csv('Sim_heat_flux_profil_t.dat', sep='\s+', header=None, names=['Time', 'z', 'Conductive', 'Advective', 'Total'])
    # Supprimer les deux premiers jours de simulation
    sim_flux = sim_flux[sim_flux['Time'] > 86400 * 2]
    sim_flux['dates'] = pd.to_datetime(sim_flux['Time'], unit='s', origin=date_simul_bg)

    # Supprimer les doublons
    sim_flux = sim_flux.drop_duplicates(subset=['Conductive', 'Advective', 'Total', 'z', 'Time'])

    # Créer une grille pour l'interpolation
    z_new = np.linspace(sim_flux['z'].min(), sim_flux['z'].max(), num=500)
    t_new = np.linspace(sim_flux['Time'].min(), sim_flux['Time'].max(), num=500)
    z_grid, t_grid = np.meshgrid(z_new, t_new)

    # Interpolation des données pour chaque type de flux
    points = np.array([sim_flux['z'], sim_flux['Time']]).T
    
    conductive_interp = griddata(points, sim_flux['Conductive'], (z_grid, t_grid), method='cubic')
    advective_interp = griddata(points, sim_flux['Advective'], (z_grid, t_grid), method='cubic')
    total_interp = griddata(points, sim_flux['Total'], (z_grid, t_grid), method='cubic')
    
    # Convertir t_new en format datetime
    t_new = pd.to_datetime(t_new, unit='s', origin=date_simul_bg)
    
    # Tracer les profils interpolés pour chaque type de flux
    fig, axes = plt.subplots(3, 1, figsize=(15, 24), sharex=True)
    
    # Conductive flux
    c1 = axes[0].contourf(t_new, z_new, conductive_interp.T, cmap='Spectral_r', levels=100)
    fig.colorbar(c1, ax=axes[0], label='Conductive flux (W/m²)')
    axes[0].set_ylabel("Depth (m)", fontsize=fontsize)
    axes[0].set_title("Profil de flux de chaleur conductif", fontsize=fontsize, pad=20)
    axes[0].grid()
    
    # Advective flux
    c2 = axes[1].contourf(t_new, z_new, advective_interp.T, cmap='Spectral_r', levels=100)
    fig.colorbar(c2, ax=axes[1], label='Advective flux (W/m²)')
    axes[1].set_ylabel("Depth (m)", fontsize=fontsize)
    axes[1].set_title("Profil de flux de chaleur advectif", fontsize=fontsize, pad=20)
    axes[1].grid()
    
    # Total flux
    c3 = axes[2].contourf(t_new, z_new, total_interp.T, cmap='Spectral_r', levels=100)
    fig.colorbar(c3, ax=axes[2], label='Total flux (W/m²)')
    axes[2].set_ylabel("Depth (m)", fontsize=fontsize)
    axes[2].set_xlabel("Date", fontsize=fontsize)
    axes[2].set_title("Profil de flux de chaleur total", fontsize=fontsize, pad=20)
    axes[2].grid()
    
    # Formater les dates sur l'axe des x
    axes[2].xaxis_date()
    fig.autofmt_xdate()

    plt.tight_layout()
    plt.show()


def plot_fluxes_timeseries(fontsize=15, date_simul_bg=None):
    """
    Plots the time series of water velocity and heat flux profiles (conductive, advective, and total) at the maximum depth (z_max).

    Parameters:
    - fontsize (int): Font size for the labels (default is 15).
    - date_simul_bg (str or datetime): Start date of the simulation in a format recognized by pandas.to_datetime.

    The function reads the water velocity and heat flux profiles from CSV files, filters the data to remove the first two days of simulation,
    and plots the time series of water velocity, conductive flux, advective flux, and total flux at the maximum depth (z_max).
    """
    # Lire les profils de vitesse
    sim_water_discharge = pd.read_csv('Sim_velocity_profil_t.dat', sep='\s+', header=None, names=['Time', 'z', 'q'])
    # delete first day of simulation
    sim_water_discharge = sim_water_discharge[sim_water_discharge['Time'] > 86400*2]
    sim_water_discharge['dates'] = pd.to_datetime(sim_water_discharge['Time'], unit='s', origin=date_simul_bg)

    # Supprimer les doublons
    sim_water_discharge = sim_water_discharge.drop_duplicates(subset=['q', 'z', 'Time'])
    # maximun value of z
    z_max=sim_water_discharge['z'].max()

    # Filtrer pour les valeurs de z supérieures ou égales à z_max
    sim_water_discharge = sim_water_discharge[sim_water_discharge['z'] >= z_max]



    # Lire les profils de flux de chaleur
    sim_flux = pd.read_csv('Sim_heat_flux_profil_t.dat', sep='\s+', header=None, names=['Time', 'z', 'Conductive', 'Advective', 'Total'])
    # Supprimer les deux premiers jours de simulation
    sim_flux = sim_flux[sim_flux['Time'] > 86400 * 2]
    sim_flux['dates'] = pd.to_datetime(sim_flux['Time'], unit='s', origin=date_simul_bg)

    # Supprimer les doublons
    sim_flux = sim_flux.drop_duplicates(subset=['Conductive', 'Advective', 'Total', 'z', 'Time'])
        # Filtrer pour les valeurs de z supérieures ou égales à z_max
    sim_flux = sim_flux[sim_flux['z'] >= z_max]

    # add 4 subplots (water velocity, conductive flux, advective flux, total flux)
    fig, axes = plt.subplots(4, 1, figsize=(10, 24), sharex=True)
    # Water velocity
    axes[0].plot(sim_water_discharge['dates'], sim_water_discharge['q'])
    axes[0].set_ylabel("Water velocity (m/s)", fontsize=fontsize)
    axes[0].set_title("Water velocity", fontsize=fontsize, pad=20)
    axes[0].grid()  
    # Conductive flux
    axes[1].plot(sim_flux['dates'], sim_flux['Conductive'])
    axes[1].set_ylabel("Conductive flux (W/m²)", fontsize=fontsize)
    axes[1].set_title("Conductive flux", fontsize=fontsize, pad=20)
    axes[1].grid()
    # Advective flux
    axes[2].plot(sim_flux['dates'], sim_flux['Advective'])
    axes[2].set_ylabel("Advective flux (W/m²)", fontsize=fontsize)
    axes[2].set_title("Advective flux", fontsize=fontsize, pad=20)
    axes[2].grid()
    # Total flux
    axes[3].plot(sim_flux['dates'], sim_flux['Total'])
    axes[3].set_ylabel("Total flux (W/m²)", fontsize=fontsize)
    axes[3].set_xlabel("Date", fontsize=fontsize)
    axes[3].set_title("Total flux", fontsize=fontsize, pad=20)
    axes[3].grid()
    # Formater les dates sur l'axe des x
    axes[3].xaxis_date()
    fig.autofmt_xdate()
    plt.tight_layout()
    plt.show()

