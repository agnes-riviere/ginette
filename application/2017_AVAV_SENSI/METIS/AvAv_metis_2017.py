#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script simple :
1. Lecture de E_geom.dat (z_haut, z_bas)
2. Lecture de E_def_mail.dat (x)
3. Garde les x uniques, associe z_haut et z_bas à chaque x
4. Trace le profil (x, z_haut) et (x, z_bas)
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import METISfunctions as mf
import Functions_mesh_metis as fm
from Functions_mesh_metis import build_gmsh_mesh_from_linestopbottom
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
geom_path = os.path.join(script_dir, '..', 'E_geom.dat')
coord_path = os.path.join(script_dir, '..', 'E_coordonnee.dat')



# 1. Lecture de E_geom.dat /home/ariviere/Programmes/ginette/application/2017_AVAV_SENSI/E_geom.dat
# trouver le fichier peu importe ou on execute le script  /home/ariviere/Programmes/ginette/application/2017_AVAV_SENSI/E_geom.dat

geom = pd.read_csv(geom_path, sep=r'\s+', header=None, names=['z_haut', 'z_bas'])




# 2. Lecture de E_def_maille.dat
# On suppose que la première colonne est x (abscisse réelle)
coordmail = pd.read_csv(coord_path, sep=r'\s+', header=None, names=['x', 'y'])

# 3. Garde les x uniques et du plus petit au plus grand
x_sorted = np.sort(np.unique(coordmail['x']))
if len(x_sorted) != len(geom):
    raise ValueError(f"Le nombre de x uniques ({len(x_sorted)}) ne correspond pas au nombre de lignes dans E_geom.dat ({len(geom)})")
z_haut_sorted = geom['z_haut'].to_numpy()
z_bas_sorted = geom['z_bas'].to_numpy()
print(f"x_sorted: {x_sorted[:5]} ... {x_sorted[-5:]}")
print(f"z_haut_sorted: {z_haut_sorted[:5]} ... {z_haut_sorted[-5:]}")
print(f"z_bas_sorted: {z_bas_sorted[:5]} ... {z_bas_sorted[-5:]}")
print(f"len(x_sorted): {len(x_sorted)}, len(z_haut_sorted): {len(z_haut_sorted)}, len(z_bas_sorted): {len(z_bas_sorted)}")

xr1=22
xr2=30
z_zns=78.
# 5. Plot
plt.figure(figsize=(10,5))
plt.plot(x_sorted, z_haut_sorted, label='z_haut (toit)', color='royalblue')
plt.plot(x_sorted, z_bas_sorted, label='z_bas (fond)', color='saddlebrown')
plt.fill_between(x_sorted, z_bas_sorted, z_haut_sorted, color='lightblue', alpha=0.3)

plt.axhline(y=z_zns, color='orange', linestyle='--', label=f'z = {z_zns} m')
plt.xlabel('x')
plt.ylabel('Altitude (m)')
plt.title('Profil extrait de E_geom.dat selon E_def_maille.dat')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Construction du maillage GMSH à partir du profil (sans cable)
mesh_filename = "AvAv_mesh.msh"
mesh_name="AvAv_mesh"

# Construction d'un rectangle à partir des 4 extrémités du profil
build_gmsh_mesh_from_linestopbottom(
    x_sorted,
    z_haut_sorted,
    z_bas_sorted[::-1],
    xr1=xr1, xr2=xr2,z_zns=z_zns,resolution_min = 0.01,resolution_max = 0.5,
    mesh_filename=mesh_filename,
    visualize=True
)
print(f"Maillage  créé et sauvegardé sous {mesh_filename}")


# gmsh to metis mesh
fm.convertMsh2Mail("AvAv_mesh")


# adapter AvAv_bck.COMM 
