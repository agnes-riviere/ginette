import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir('/home/ariviere/Programmes/ginette/application/1D_col')

# Charger les données depuis le fichier "Sim_pressure_profil_t.dat"
data = np.loadtxt("Sim_pressure_profil_t.dat")
data_sat = np.loadtxt("S_saturation_profil_t.dat")
data_perm = np.loadtxt("S_permeabilite_t.dat")
data_temp = np.loadtxt("Sim_temperature_profil_t.dat")

# Extraire les valeurs de z, t et pression depuis les colonnes du tableau "data"
z = data[:, 1]
t = data[:, 0]
pression = data[:, 3]
# Extraire les valeurs de z, t et saturation depuis les colonnes du tableau "data_sat"
z_sat = data_sat[:, 1]
t_sat = data_sat[:, 0]
sat = data_sat[:, 2]
# Extraire les valeurs de z, t et saturation depuis les colonnes du tableau "data_perm"
z_perm = data_perm[:, 1]
t_perm = data_perm[:, 0]
perm = data_perm[:, 2]
# Extraire les valeurs de z, t et température depuis les colonnes du tableau "data_temp"
z_temp = data_temp[:, 1]
t_temp = data_temp[:, 0]
temperature = data_temp[:, 2]

# Create subplots for pressure, saturation, permeability, and temperature data
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

# Plot pressure data
sc1 = ax1.scatter(t, z, c=pression, cmap='viridis', label='H (m)')
ax1.set_xlabel('Temps (s)')
ax1.set_ylabel('Profondeur (m)')
ax1.set_title('H en fonction de z et de t')
cbar1 = plt.colorbar(sc1, ax=ax1)
cbar1.set_label('H (m)')

# Plot saturation data
sc2 = ax2.scatter(t_sat, z_sat, c=sat, cmap='viridis', label='sat')
ax2.set_xlabel('Temps (s)')
ax2.set_ylabel('Profondeur (m)')
ax2.set_title('Saturation en fonction de z et de t')
cbar2 = plt.colorbar(sc2, ax=ax2)
cbar2.set_label('Saturation')

# Apply logarithmic scaling to permeability values
perm_log = np.log10(perm)

# Plot permeability data on a logarithmic scale
sc3 = ax3.scatter(t_perm, z_perm, c=perm_log, cmap='viridis', label='perm')
ax3.set_xlabel('Temps (s)')
ax3.set_ylabel('Profondeur (m)')
ax3.set_title('Perméabilité en fonction de z et de t (échelle log)')
cbar3 = plt.colorbar(sc3, ax=ax3)
cbar3.set_label('log(Perméabilité)')

# Plot temperature data
sc4 = ax4.scatter(t_temp, z_temp, c=temperature, cmap='inferno', label='temp')
ax4.set_xlabel('Temps (s)')
ax4.set_ylabel('Profondeur (m)')
ax4.set_title('Température en fonction de z et de t')
cbar4 = plt.colorbar(sc4, ax=ax4)
cbar4.set_label('Température (°C)')

plt.tight_layout()
plt.show()
