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
import re
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

# adapte les points des contours verticaux dans le .COMM au maillage courant
comm_hydro_path = os.path.join(script_dir, "AvAv_hydro.COMM")
if os.path.isfile(comm_hydro_path):
    x_min = float(np.min(x_sorted))
    x_max = float(np.max(x_sorted))
    z_top_left = float(z_haut_sorted[np.argmin(x_sorted)])
    z_top_right = float(z_haut_sorted[np.argmax(x_sorted)])
    x_mid = 0.5 * (x_min + x_max)
    z_top_mid = float(np.interp(x_mid, x_sorted, z_haut_sorted))
    z_bot_left = float(z_bas_sorted[np.argmin(x_sorted)])
    z_bot_right = float(z_bas_sorted[np.argmax(x_sorted)])
    z_mid_left = 0.5 * (z_top_left + z_bot_left)
    z_mid_right = 0.5 * (z_top_right + z_bot_right)

    with open(comm_hydro_path, "r", encoding="utf-8") as f_comm:
        comm_txt = f_comm.read()

    # côté x = xmin (ligne verticale) : remplace toute ligne "points ... #coté droit vertical"
    comm_txt = re.sub(
        r"points\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+#coté droit vertical",
        f"points  {x_min:.3f} {z_top_left:.3f}  {x_min:.3f} {z_bot_left:.3f}  {x_min:.3f} {z_mid_left:.3f} #coté droit vertical",
        comm_txt,
    )

    # côté x = xmax (ligne verticale) : remplace toute ligne "points ... #coté gauche vertical"
    comm_txt = re.sub(
        r"points\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+#coté gauche vertical",
        f"points {x_max:.3f} {z_bot_right:.3f}  {x_max:.3f} {z_top_right:.3f}  {x_max:.3f} {z_mid_right:.3f} #coté gauche vertical",
        comm_txt,
    )

    # surface topo: remplace toute ligne "points ... # surface topo"
    comm_txt = re.sub(
        r"points\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+# surface topo",
        f"points {x_max:.3f} {z_top_right:.3f}  {x_min:.3f} {z_top_left:.3f}  {x_mid:.3f} {z_top_mid:.3f} # surface topo",
        comm_txt,
    )

    with open(comm_hydro_path, "w", encoding="utf-8") as f_comm:
        f_comm.write(comm_txt)

    print(
        "AvAv_hydro.COMM adapté au maillage pour les contours verticaux: "
        f"xmin={x_min:.3f}, xmax={x_max:.3f}, "
        f"zbot~{min(z_bot_left, z_bot_right):.3f}, ztop~{max(z_top_left, z_top_right):.3f}"
    )

# adapter AvAv_bck.COMM 
# utilise : /home/ariviere/Programmes/ginette/application/2017_AVAV_SENSI/E_Id_river.dat et /home/ariviere/Programmes/ginette/application/2017_AVAV_SENSI/E_coordonnee.dat
# donne moi les coordonnées de la rivière (x_riv, z_riv) à partir de E_Id_river.dat et E_coordonnee.dat
riv_id_path = os.path.join(script_dir, '..', 'E_Id_river.dat')
riv_id = pd.read_csv(riv_id_path, sep=r'\s+', header=None, names=['id'])
coord_path = os.path.join(script_dir, '..', 'E_coordonnee.dat')
coord = pd.read_csv(coord_path, sep=r'\s+', header=None, names=['x', 'z'])
# numéro de ligne coord = id maille retrouver les coordonnées des riv_id
# E_Id_river.dat provient de Fortran (indexation 1-based): on convertit en index pandas (0-based)
river_ids = riv_id['id'].to_numpy(dtype=int)
if river_ids.min() < 1 or river_ids.max() > len(coord):
    raise ValueError(
        f"IDs rivière hors bornes: min={river_ids.min()}, max={river_ids.max()}, "
        f"nb_coord={len(coord)}"
    )

river_index = river_ids - 1

riv_coords = coord.iloc[river_index].reset_index(drop=True)
x_riv = riv_coords['x'].to_numpy()
z_riv = riv_coords['z'].to_numpy()
print(f"Coordonnées de la rivière (x_riv, z_riv) : {list(zip(x_riv, z_riv))[:5]} ... {list(zip(x_riv, z_riv))[-5:]}")


# fill /home/ariviere/Programmes/ginette/application/2017_AVAV_SENSI/METIS/trace_riv_berge avec les coordonnées de gmsh pour les  points de la rivières c'est à dire entre les x_min et x_max de la rivière et les noeuds de la surface du maillage compris entre les deux x_min et x_max de la rivière
mail_path = os.path.join(script_dir, "AvAv_mesh.mail")
trace_path = os.path.join(script_dir, "trace_riv_berge")

with open(mail_path, "r", encoding="utf-8") as f_mail:
    header = f_mail.readline().split()
    n_nodes = int(header[0])
    node_lines = [f_mail.readline() for _ in range(n_nodes)]

node_coords = np.array([
    [float(parts[0]), float(parts[1])]
    for parts in (line.split() for line in node_lines)
], dtype=float)

node_x = node_coords[:, 0]
node_z = node_coords[:, 1]
node_ids = np.arange(1, n_nodes + 1, dtype=int)  # base 1 METIS

x_min_riv = float(np.min(x_riv))
x_max_riv = float(np.max(x_riv))

# altitude de surface théorique au x des noeuds du maillage
z_top_interp = np.interp(node_x, x_sorted, z_haut_sorted)

# tolérance verticale sur la détection de la surface
dz_tol = 5e-3
mask_surface_river = (
    (node_x >= x_min_riv) &
    (node_x <= x_max_riv) &
    (np.abs(node_z - z_top_interp) <= dz_tol)
)

surface_river_nodes = node_ids[mask_surface_river]
surface_river_x = node_x[mask_surface_river]

if surface_river_nodes.size == 0:
    raise ValueError(
        "Aucun noeud de surface trouvé pour trace_riv_berge "
        f"(x in [{x_min_riv:.3f}, {x_max_riv:.3f}], dz_tol={dz_tol})"
    )

# ordre gauche -> droite
order = np.argsort(surface_river_x)
surface_river_nodes = surface_river_nodes[order]

# écriture au format attendu par METIS include: lignes préfixées par 'no'
with open(trace_path, "w", encoding="utf-8") as f_trace:
    chunk_size = 10
    for i in range(0, len(surface_river_nodes), chunk_size):
        chunk = surface_river_nodes[i:i + chunk_size]
        f_trace.write("   no  " + " ".join(str(n) for n in chunk) + "\n")

print(
    f"trace_riv_berge généré: {trace_path} | "
    f"{len(surface_river_nodes)} noeuds de surface entre x={x_min_riv:.3f} et x={x_max_riv:.3f}"
)

# figure verification
plt.figure(figsize=(10, 5))
plt.plot(node_x, node_z, 'o', markersize=2, label='Noeuds du maillage', alpha=0.5)
plt.plot(x_sorted, z_haut_sorted, label='z_haut (toit)', color='royalblue')
plt.plot(x_sorted, z_bas_sorted, label='z_bas (fond)', color='saddlebrown')
plt.fill_between(x_sorted, z_bas_sorted, z_haut_sorted, color='lightblue', alpha=0.3)
plt.scatter(surface_river_x, node_z[mask_surface_river], color='red', label='Noeuds de surface rivière', zorder=5)
plt.axhline(y=z_zns, color='orange', linestyle='--', label=f'z = {z_zns} m')
plt.xlabel('x')
plt.ylabel('Altitude (m)')
plt.title('Vérification des noeuds de surface pour trace_riv_berge')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

input_list = ['E_chargeT_RD.dat','E_chargeT_RG.dat','E_chargeT_Riv.dat','E_tempT_RD.dat','E_tempT_RG.dat','E_tempT_Riv.dat']
output_list = ['Ch_RD.txt','Ch_RG.txt','Ch_R.txt','temp_RD.txt','temp_RG.txt','temp_R.txt']
mf.creation_fichier_txt(input_list=input_list,output_list=output_list)

# conversion des conditions limites en .fic en utilisan t les fonctions de METISfunctions
boundary_files_h = fm.run_cree_fic3(
    input_r_txt="Ch_R.txt",
    output_r_fic="Ch_R.fic",
    input_rd_txt="Ch_RD.txt",
    output_rd_fic="Ch_RD.fic",
    output_rg_fic="Ch_RG.fic"
)

# fichier attendu par AvAv_bck.COMM : temperature contour fichier = T_R.fic
boundary_files_t = fm.run_cree_fic3(
    input_r_txt="temp_R.txt",
    output_r_fic="T_R.fic",
    input_rd_txt="temp_RD.txt",
    output_rd_fic="T_RD.fic",
    output_rg_fic="T_RG.fic"
)

print(f"Fichiers de conditions limites hydrauliques créés : {boundary_files_h}")
print(f"Fichiers de conditions limites thermiques créés : {boundary_files_t}")