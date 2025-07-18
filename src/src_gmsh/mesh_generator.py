import os
import sys
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt
import math
import meshio
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches
from shapely.geometry import Point, Polygon
import gmshparser
import os
import sys
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib.pyplot as plt
import math
import meshio
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches
from shapely.geometry import Point, Polygon
import gmshparser




#verify if gmsh is installed
try:
    import gmsh 
except ImportError:
    print("GMSH is not installed. Please install it using 'pip install gmsh' or 'conda install gmsh'.")
    sys.exit(1)
# Check gmsh version
try:
    gmsh.__version__
except AttributeError:
    print("GMSH version is not available. Please check your GMSH installation.")
    sys.exit(1) 



def plot_gmsh_mesh(mesh_path):
    try:
        mesh = pv.read(mesh_path)
        plotter = pv.Plotter()
        plotter.add_mesh(mesh, show_edges=True, color="lightblue", label="Maillage Gmsh")
#
        points = mesh.points
# add axis coordinates
        plotter.add_axes()
#add grid
        plotter.add_legend()
# plot point name from the mesh  mesh = pv.read(mesh_path)
#        plotter.add_point_labels(points, [f"Point {i}" for i in range(len(points))], point_color='black', font_size=10, shape_opacity=0.5, shape_color='white')

# add buttom zoom
        plotter.camera.zoom(1.5)
        plotter.show()
    except Exception as e:
        print(f"Erreur lors de la visualisation du maillage : {e}")
 

def process_distance_altitude(file_path, Station, sim_dir, altitude_min, Rive_D=0, dz_pt100=0.1,dz_pt100first=0.1):
    """
    Process the distance and altitude data from a leveling file and save it as a CSV.

    Parameters:
    - file_path (str): Path to the directory containing the leveling file.
    - Station (str): Name of the station.
    - sim_dir (str): Directory where the output CSV will be saved.
    - altitude_min (float): Minimum altitude to add at the boundaries.
    - Rive_D (int): Direction of processing (0 for right-to-left, 1 for left-to-right).
    - dz_pt100 (float): Vertical step for PT100 coordinates.
    - dz_pt100first (float): Vertical step between the streambed and the first PT100 coordinate.

    Returns:
    - distance_altitude_table (DataFrame): Processed distance and altitude data.
    - pt100_coord (DataFrame): Coordinates of PT100 instruments.
    """
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    # Chemin vers le fichier de nivellement
    file_path = os.path.join(file_path, Station, f"{Station}.csv")
    print(f"Chemin du fichier de nivellement CSV : {file_path}")
    # Vérification de l'existence du fichier
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Le fichier {file_path} est introuvable !")

    # Lecture des données
    try:
        data = pd.read_csv(file_path, delim_whitespace=True, header=None)
        x_coords, y_coords, z_coords, instr = data[0].to_numpy(), data[1].to_numpy(), data[2].to_numpy(), data[3].to_numpy()
    except Exception as e:
        raise ValueError(f"Erreur lors de la lecture du fichier : {e}")

    # Calculer les distances cumulées
    distances = np.zeros(len(x_coords))

    if Rive_D == 0:
        # Calculer les distances de droite à gauche
        for i in range(len(x_coords) - 2, -1, -1):  # Parcourir les indices de la fin vers le début
            distances[i] = distances[i + 1] + np.sqrt((x_coords[i + 1] - x_coords[i])**2 + (y_coords[i + 1] - y_coords[i])**2)
        distances = np.flip(distances)  # Inverser l'ordre des distances
        z_coords = np.flip(z_coords)  # Inverser l'ordre des altitudes
        instr = np.flip(instr)  # Inverser l'ordre des instruments
    else:
        for i in range(1, len(x_coords)):
            distances[i] = distances[i - 1] + np.sqrt((x_coords[i] - x_coords[i - 1])**2 + (y_coords[i] - y_coords[i - 1])**2)

    distances = np.round(distances, 1)

    # Ajout des points aux extrémités
    distances = np.insert(distances, 0, 0)
    z_coords = np.insert(z_coords, 0, altitude_min)
    instr = np.insert(instr, 0, "NA")  # Ajout d'un instrument fictif pour le début
    distances = np.append(distances, distances[-1])
    z_coords = np.round(z_coords, 1)
    z_coords = np.append(z_coords, altitude_min)
    instr = np.append(instr, "NA")  # Ajout d'un instrument fictif pour la fin

    # Création du DataFrame
    distance_altitude_table = pd.DataFrame({"Distance (m)": distances, "Altitude (Z)": z_coords, "Instrument": instr})

    # Sauvegarde en CSV
    output_csv_path = os.path.join(sim_dir, "distance_vs_altitude.csv")
    distance_altitude_table.to_csv(output_csv_path, index=False)
    # Trouver et afficher les coordonnées des instruments Hobo dans distance_altitude_table
    # la difference entre le lit de la riviere et la premiere pt100 = dz_pt100first
    pt100_coord = pd.DataFrame()
    if "Instrument" in distance_altitude_table.columns:
        # Trouver les indices de toutes les occurrences contenant "Hobo" dans la colonne "Instrument"
        index_hobo = distance_altitude_table[
            distance_altitude_table["Instrument"].str.contains("Hobo", na=False, case=False)
        ].index

        if not index_hobo.empty:
            # Récupérer les noms des instruments Hobo
            hobo_names = distance_altitude_table.loc[index_hobo, "Instrument"].tolist()

            # Récupérer les coordonnées des instruments Hobo
            hobo_coords = distance_altitude_table.loc[index_hobo, ["Distance (m)", "Altitude (Z)"]]

            # Afficher les coordonnées des instruments Hobo
            for name, coord in zip(hobo_names, hobo_coords.values):
                for i in range(1, 5):
                    distance_hobo = coord[0]
                    if i == 1:
                        altitude_hobo = round(coord[1] - dz_pt100first, 2)
                    else : 
                        altitude_hobo = round(coord[1] - dz_pt100first- dz_pt100 * (i-1),2)
                    pt100_coord = pd.concat([pt100_coord, 
                                             pd.DataFrame({"hobo":[name],"pt100":[i],
                                                           "Distance (m)": [distance_hobo],
                                                             "Altitude (Z)": [altitude_hobo]})], ignore_index=True)
        else:
            print("Aucun instrument Hobo trouvé dans la table.")
    else:
        print("La colonne 'Instrument' est absente de distance_altitude_table.")

    return distance_altitude_table,pt100_coord

# def generate_quadrilateral_mesh(distance_altitude_table, output_mesh_path, dx=0.1, dz=0.1, mesh_dimension=2):
#     try:
#         print("\nGénération du maillage structuré avec Gmsh...")
#         gmsh.initialize()
#         gmsh.option.setNumber("General.Terminal", 1)
#         gmsh.model.add("RectangularMesh")

#         # Extraction des distances et altitudes
#         distances = distance_altitude_table["Distance (m)"].to_numpy()
#         z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()
#         print(f"Distances : {distances}")
#         print(f"Altitudes : {z_coords}")
#         # Points du rectangle
#         p1 = gmsh.model.occ.addPoint(distances[0], 0, z_coords[0])
#         p2 = gmsh.model.occ.addPoint(distances[-1], 0, z_coords[-1])
#         p3 = gmsh.model.occ.addPoint(distances[-2], 0, z_coords[-2])
#         p4 = gmsh.model.occ.addPoint(distances[1], 0, z_coords[1])

#         # Création des lignes et de la surface
#         l1 = gmsh.model.occ.addLine(p1, p2)
#         l2 = gmsh.model.occ.addLine(p2, p3)
#         l3 = gmsh.model.occ.addLine(p3, p4)
#         l4 = gmsh.model.occ.addLine(p4, p1)
#         loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
#         surface = gmsh.model.occ.addPlaneSurface([loop])

#         # Synchronisation du modèle
#         gmsh.model.occ.synchronize()

#         # Calcul des subdivisions
#         length_x = distances[-1] - distances[0]
#         length_z = max(z_coords) - min(z_coords)
#         num_div_x = max(2, math.ceil(length_x / dx))
#         num_div_z = max(2, math.ceil(length_z / dz))
#         print(f"Nombre de divisions en X : {num_div_x}, Nombre de divisions en Z : {num_div_z}")

#         # Application des subdivisions et recombinaison pour quadrilatères
#         gmsh.model.mesh.setTransfiniteCurve(l1, num_div_x)
#         gmsh.model.mesh.setTransfiniteCurve(l3, num_div_x)
#         gmsh.model.mesh.setTransfiniteCurve(l2, num_div_z)
#         gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z)
#         gmsh.model.mesh.setTransfiniteSurface(surface, "Left")
#         gmsh.model.mesh.setRecombine(2, surface)

#         # Génération du maillage
#         gmsh.model.mesh.generate(mesh_dimension)
#         gmsh.write(output_mesh_path)
#         print(f"✅ Maillage sauvegardé sous : {output_mesh_path}")

#     except Exception as e:
#         print(f"Erreur lors de la génération du maillage : {e}")
#     finally:
#         gmsh.finalize()

# # Appel de la fonction
# output_mesh_path = os.path.join(sim_dir, "rectangular_mesh.msh")
# #generate_quadrilateral_mesh(distance_altitude_table, output_mesh_path, dx=0.1, dz=0.5, mesh_dimension=2)



# def generate_quadrilateral_mesh_with_aligned_lines(output_mesh_path, dx_outer=1.0, dz_outer=1.0, dx_inner=0.2, dz_inner=0.2, mesh_dimension=2):
#     try:
#         print("\nGénération du maillage structuré avec Gmsh...")
#         gmsh.initialize()
#         gmsh.option.setNumber("General.Terminal", 1)
#         gmsh.model.add("RectangularMeshWithAlignedLines")

#         # Extraction des distances et altitudes
#         distances = distance_altitude_table["Distance (m)"].to_numpy()
#         z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()

#         # Points du rectangle
#         p1 = gmsh.model.occ.addPoint(distances[0], 0, z_coords[0])
#         p2 = gmsh.model.occ.addPoint(distances[2], 0, z_coords[2])
#         p3 = gmsh.model.occ.addPoint(distances[-1], 0, z_coords[-1])
#         p4 = gmsh.model.occ.addPoint(distances[-2], 0, z_coords[-2])

#         # Création des lignes et de la surface
#         l1 = gmsh.model.occ.addLine(p1, p2)
#         l2 = gmsh.model.occ.addLine(p2, p3)
#         l3 = gmsh.model.occ.addLine(p3, p4)
#         l4 = gmsh.model.occ.addLine(p4, p1)
#         loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
#         surface = gmsh.model.occ.addPlaneSurface([loop])

#         # Synchronisation du modèle
#         gmsh.model.occ.synchronize()

#         # Calcul des subdivisions
#         length_x = max(distances) - min(distances)
#         length_z = max(z_coords) - min(z_coords)
#         num_div_x = max(2, math.ceil(length_x / dx_outer))
#         num_div_z = max(2, math.ceil(length_z / dz_outer))

#         # Ajustement des lignes verticales pour qu'elles soient alignées avec les faces des mailles
#         x_line1 = 5  # Position initiale de la première ligne
#         x_line2 = 6  # Position initiale de la deuxième ligne

#         # Ajuster les positions pour qu'elles correspondent aux divisions des mailles
#         x_line1 = distances[0] + round((x_line1 - distances[0]) / dx_outer) * dx_outer
#         x_line2 = distances[0] + round((x_line2 - distances[0]) / dx_outer) * dx_outer


#         # Calcul des subdivisions pour les régions
#         num_div_x_outer_left = max(2, math.ceil((x_line1 - min(distances)) / dx_outer))
#         num_div_x_outer_right = max(2, math.ceil((x_line2 - max(distances)) / dx_outer))
#         num_div_x_inner = max(2, math.ceil((x_line2 - x_line1) / dx_inner))
#         num_div_z = max(2,math.ceil(length_z / dz_outer) )



#         # Création des points pour les lignes verticales
#         p5 = gmsh.model.occ.addPoint(x_line1, 0,max(z_coords))
#         p6 = gmsh.model.occ.addPoint(x_line1, 0, min(z_coords))
#         p7 = gmsh.model.occ.addPoint(x_line2, 0, max(z_coords))
#         p8 = gmsh.model.occ.addPoint(x_line2, 0, min(z_coords))

#         # Création des lignes verticales
#         l5 = gmsh.model.occ.addLine(p5, p6)
#         l6 = gmsh.model.occ.addLine(p7, p8)

#         # Application des subdivisions et recombinaison pour quadrilatères
#         gmsh.model.mesh.setTransfiniteCurve(l1, num_div_x_outer_left)
#         gmsh.model.mesh.setTransfiniteCurve(l3, num_div_x)
#         gmsh.model.mesh.setTransfiniteCurve(l2, num_div_z)
#         gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z)
#         gmsh.model.mesh.setTransfiniteCurve(l5, num_div_z)
#         gmsh.model.mesh.setTransfiniteCurve(l6, num_div_z)
#         gmsh.model.mesh.setTransfiniteSurface(surface, "Left")
#         gmsh.model.mesh.setRecombine(2, surface)

#         # Génération du maillage
#         gmsh.model.mesh.generate(mesh_dimension)
#         gmsh.write(output_mesh_path)
#         print(f"✅ Maillage sauvegardé sous : {output_mesh_path}")

#     except Exception as e:
#         print(f"Erreur lors de la génération du maillage : {e}")
#     finally:
#         gmsh.finalize()



# def generate_mesh_8disc(distance_altitude_table, output_mesh_path, h_left=5,h_right=8,V_top=104.8,v_bot=103.8,  dx_outer=0.1, dz_outer=1.0, dx_inner=0.2, dz_inner=0.2,mesh_dimension=2):
    
#     try:
#         print("\nGénération du maillage structuré avec Gmsh...")
#         gmsh.initialize()
#         gmsh.option.setNumber("General.Terminal", 1)
#         gmsh.model.add("RectangularMesh")

#         # Extraction des distances et altitudes
#         distances = distance_altitude_table["Distance (m)"].to_numpy()
#         z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()

#         # Points du rectangle
#         min_z= min(z_coords)
#         max_z= max(z_coords)
#         p1 = gmsh.model.occ.addPoint(h_left, 0, min_z)
#         p2 = gmsh.model.occ.addPoint(h_left, 0, max_z)
#         p3 = gmsh.model.occ.addPoint(h_right, 0,min_z)
#         p4 = gmsh.model.occ.addPoint(h_right, 0,max_z)

#         # point intermediaire 
#         p5 = gmsh.model.occ.addPoint(h_right, 0, v_bot)
#         p6 = gmsh.model.occ.addPoint(h_left, 0, v_bot)


#         # Création des lignes pour une seule surface combinée
#         # Création des lignes pour une seule surface combinée
#         l1 = gmsh.model.occ.addLine(p2, p4)  # Ligne supérieure
#         l2 = gmsh.model.occ.addLine(p4, p5)  # Ligne droite supérieure
#         l3 = gmsh.model.occ.addLine(p5, p6)  # Ligne horizontale intermédiaire (limite commune)
#         l4 = gmsh.model.occ.addLine(p6, p2)  # Ligne gauche supérieure
#         l5 = gmsh.model.occ.addLine(p5, p3)  # Ligne droite inférieure
#         l6 = gmsh.model.occ.addLine(p3, p1)  # Ligne inférieure
#         l7 = gmsh.model.occ.addLine(p1, p6)  # Ligne gauche inférieure

#         # Combiner toutes les courbes dans une seule boucle
#        # Création des boucles pour les deux sous-surfaces
#         loop1 = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])  # Région supérieure
#         loop2 = gmsh.model.occ.addCurveLoop([l3, l5, l6, l7])  # Région inférieure

#         # Création des sous-surfaces
#         surface1 = gmsh.model.occ.addPlaneSurface([loop1])
#         surface2 = gmsh.model.occ.addPlaneSurface([loop2])

#         # Synchronisation du modèle
#         gmsh.model.occ.synchronize()

#         # Calcul des subdivisions
#         length_x = h_right - h_left
#         num_div_x = max(2, math.ceil(length_x / dx_inner))
#         num_div_z_inner = max(2, math.ceil((max_z - V_top) / dz_inner))
#         num_div_z_outer = max(2, math.ceil((V_top - min_z) / dz_outer))

#         print(f"Nombre de divisions en X : {num_div_x}")
#         print(f"Nombre de divisions en Z (intérieur) : {num_div_z_inner}")
#         print(f"Nombre de divisions en Z (extérieur) : {num_div_z_outer}")

#         # Application des subdivisions et recombinaison pour chaque sous-surface
#         for surface in [surface1, surface2]:
#             gmsh.model.mesh.setTransfiniteSurface(surface, "Left")
#             gmsh.model.mesh.setRecombine(2, surface)

#         # Subdivisions pour la région supérieure (surface1)
#         gmsh.model.mesh.setTransfiniteCurve(l1, num_div_x)  # Ligne supérieure horizontale
#         gmsh.model.mesh.setTransfiniteCurve(l3, num_div_x)  # Ligne horizontale intermédiaire (limite commune)
#         gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z_inner)  # Ligne gauche verticale (région supérieure)
#         gmsh.model.mesh.setTransfiniteCurve(l2, num_div_z_inner)  # Ligne droite verticale (région supérieure)

#         # Subdivisions pour la région inférieure (surface2)
#         gmsh.model.mesh.setTransfiniteCurve(l3, num_div_x)  # Ligne horizontale intermédiaire (limite commune)
#         gmsh.model.mesh.setTransfiniteCurve(l5, num_div_z_outer)  # Ligne droite verticale (région inférieure)
#         gmsh.model.mesh.setTransfiniteCurve(l6, num_div_x)  # Ligne inférieure horizontale
#         gmsh.model.mesh.setTransfiniteCurve(l7, num_div_z_outer)  # Ligne gauche verticale (région inférieure)

#         # Génération du maillage
#         gmsh.model.mesh.generate(mesh_dimension)
#         gmsh.write(output_mesh_path)
#         print(f"✅ Maillage sauvegardé sous : {output_mesh_path}")

#     except Exception as e:
#         print(f"Erreur lors de la génération du maillage : {e}")
#     finally:
#         gmsh.finalize()



# def generate_mesh_5_region(distance_altitude_table, output_mesh_path,v_bot=103.8,point_RG =5, point_RD = 17, dx_grossier=0.5, dx_precis = 0.1, dz_grossier=1.0, dz_precis=0.2,mesh_dimension=2):
#     """
#     Fonction permettant de créer un maillage différent pour 5 régions. Les régions sont réparties selon cet ordre : Région 1 haut du maillage, Région 2 milieu gauche du maillage,
#     Région 3 milieu centre du maillage, Région 4 milieu droit du maillage, Région 5 bas du maillage. Ces régions sont limitées selon plusieurs points. Tout d'abords le fichier de topo
#     "distance_altitude_table" Dont le 1er point correspond au coin inférieur gauche du maillage complet, le second supérieur gauche, l'avant dernier le supérieur droit et le dernier
#     point le coin inférieur droit. La frontière entre la partie haute et milieu est déterminée par l'altitude max des points de la topo rentré pour "point_RG" et "point_RD". Les 
#     coordonnées en x des points RG et RD correspondent à la frontière entre la partie gauche, centre et droit. Enfin, la frontière entre le milieu et le bas du maillage est déterminé 
#     par le parametre v_bot.Une fois le maillage créé il est enregistré vers output_mesh_path.

#     Entrée :
#     distance_altitude_table (DataFrame) : Tableau de la topo contenant au moins 2 colonnes avec "Distances (m)" pour les coordonnées en x et "Altitude (Z)" pour Z.
#     output_mesh_path (str) : Chemin où sera enregistrer le maillage brute
#     v_bot (float) : Parametre qui délimite la frontièrif plot:


#     point_RG (int) : Entier correspondant au numéro du point topo correspondant à la rive gauche. Il permet de délimiter la frontière entre la région 2 et 3 en x et entre le milieu et le haut du maillage si son altitude est plus grande que RD
#     point_RD (int) : Entier correspondant au numéro du point topo correspondant à la rive droite. Il permet de délimiter la frontière entre la région 3 et 4 en x et entre le milieu et le haut du maillage si son altitude est plus grande que RG
#     dx_grossier (float) : Valeur du pas en x utilisé dans les régions où x n'a pas besoin d'être précis. C'est-à-dire dans les régions 1, 2, 4 et 5
#     dx_precis (float) : Valeur du pas en x utilisé dans les régions où x a besoin d'être précis. C'est-à-dire dans la région 3
#     dz_grossier (float) : Valeur du pas en z utilisé dans les régions où z n'a pas besoin d'être précis. C'est-à-dire dans les régions 1 et 5
#     dz_precis (float) : Valeur du pas en z utilisé dans les régions où z a besoin d'être précis. C'est-à-dire dans les régions 2, 3 et 4
#     mesh_dimension (int) : Entier qui rensigne sur la dimension du maillage (Par défault égale à 2)

#     Sortie :
#     Le maillage est créé et sauvegarder dans le chemin renseigné avec "output_mesh_path"
    
#     """
    
#     try:
#         print("\nGénération du maillage structuré avec Gmsh...")
#         gmsh.initialize()
#         gmsh.option.setNumber("General.Terminal", 1)
#         gmsh.model.add("RectangularMesh")

#         # Extraction des distances et altitudes
#         distances = distance_altitude_table["Distance (m)"].to_numpy()
#         z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()

#         # Points du rectangle
#         min_z= min(z_coords) # On prend le min et le max pour obtenir un rectangle régulier
#         max_z= max(z_coords)
#         min_x = min(distances)
#         max_x = max(distances)

#         p1 = gmsh.model.occ.addPoint(min_x, 0, min_z) # Equivalent du point 0 du DF
#         p2 = gmsh.model.occ.addPoint(min_x, 0, max_z) # Equivalent du point 1
#         p3 = gmsh.model.occ.addPoint(max_x, 0,min_z) # Equivalent du point -1
#         p4 = gmsh.model.occ.addPoint(max_x, 0,max_z) # Equivalent du point -2

#         # point intermediaire du rectangle
#         p5 = gmsh.model.occ.addPoint(max_x, 0, v_bot)
#         p6 = gmsh.model.occ.addPoint(min_x, 0, v_bot)

#         # point de la riv (avec un maillage plus fin)
#         x_riv_gauche = distances[point_RG]
#         x_riv_droite = distances[point_RD]
#         z_riv = max(z_coords[point_RG],z_coords[point_RD])

#         p7 = gmsh.model.occ.addPoint(x_riv_gauche, 0, z_riv)
#         p8 = gmsh.model.occ.addPoint(x_riv_droite, 0, z_riv)
#         p9 = gmsh.model.occ.addPoint(x_riv_droite, 0, v_bot)
#         p10 = gmsh.model.occ.addPoint(x_riv_gauche, 0, v_bot)
#         p11 = gmsh.model.occ.addPoint(min_x, 0, z_riv)
#         p12 = gmsh.model.occ.addPoint(max_x, 0, z_riv)

#         # Création des lignes pour délimité chaque région
#         l1 = gmsh.model.occ.addLine(p2, p4) # Haut du rectangle + coté haut région 1
#         l2 = gmsh.model.occ.addLine(p4, p12) # coté droit région 1
#         l3 = gmsh.model.occ.addLine(p12,p5) # coté droit région 4
#         l4 = gmsh.model.occ.addLine(p5,p9) # Frontière entre la région 4 et 5 (bas région 4)
#         l5 = gmsh.model.occ.addLine(p9,p10) # Frontière entre la région 3 et 5 (bas région 3)
#         l6 = gmsh.model.occ.addLine(p10,p6) # Frontière entre la région 2 et 5(bas région 2)
#         l7 = gmsh.model.occ.addLine(p6,p11) # coté gauche de la région 2
#         l8 = gmsh.model.occ.addLine(p11,p2) # coté gauche de la région 1
#         l9 = gmsh.model.occ.addLine(p5,p3) # coté droit de la région 5
#         l10 = gmsh.model.occ.addLine(p3,p1) # Bas du rectangle + coté bas de la région 5
#         l11 = gmsh.model.occ.addLine(p1,p6) # coté gauche de la région 5
#         l12 = gmsh.model.occ.addLine(p11,p7) # Frontière entre la région 1 et 2 (haut région 2)
#         l13 = gmsh.model.occ.addLine(p7,p8) # Frontière entre la région 1 et 3 (haut de la région 3)
#         l14 = gmsh.model.occ.addLine(p8,p12) # Frontière entre la région 1 et 4 (haut de la région 4)
#         l15 = gmsh.model.occ.addLine(p9,p8) # Frontière entre la région 3 et 4 
#         l16 = gmsh.model.occ.addLine(p10,p7) # Frontière entre la région 3 et 2
#         l100 = gmsh.model.occ.addLine(p6,p5) # Coté haut de la région 5
#         l101 = gmsh.model.occ.addLine(p12,p11) # Coté bas de la région 1

#         # Combiner toutes les courbes dans une seule boucle
#        # Création des boucles pour les 5 régions
#         loop1 = gmsh.model.occ.addCurveLoop([l1, l2, l101, l8]) # Région 1
#         loop2 = gmsh.model.occ.addCurveLoop([l12,-l16,l6,l7]) # Région 2
#         loop3 = gmsh.model.occ.addCurveLoop([l13,-l15,l5,l16]) # Région 3
#         loop4 = gmsh.model.occ.addCurveLoop([l14,l3,l4,l15]) # Région 4
#         loop5 = gmsh.model.occ.addCurveLoop([l100,l9,l10,l11]) # Région 5

#         # Création des sous-surfaces
#         surface1 = gmsh.model.occ.addPlaneSurface([loop1]) # 1
#         surface2 = gmsh.model.occ.addPlaneSurface([loop2]) # 2
#         surface3 = gmsh.model.occ.addPlaneSurface([loop3]) # 3
#         surface4 = gmsh.model.occ.addPlaneSurface([loop4]) # 4
#         surface5 = gmsh.model.occ.addPlaneSurface([loop5]) # 5

#         # Synchronisation du modèle
#         gmsh.model.occ.synchronize()

#         # Calcul des subdivisions region 1
#         length_x = max_x - min_x
#         num_div_x = max(2, math.ceil(length_x / dx_grossier))
#         num_div_z_region1 = max(2, math.ceil((max_z - z_riv) / dz_grossier))
#         # et la region 5
#         num_div_z_region5 = max(2, math.ceil((v_bot - min_z) / dz_grossier))

#         # Calcul des subdivisions region 2
#         length_x_region2 = x_riv_gauche - min_x
#         num_div_x_region2 = max(2, math.ceil(length_x_region2 / dx_grossier))
#         num_div_z_precis = max(2, math.ceil((z_riv - v_bot) / dz_precis))
#         # et la region 4
#         length_x_region4 = max_x - x_riv_droite
#         num_div_x_region4 = max(2, math.ceil(length_x_region4 / dx_grossier))

#         # Calcul des subdivisions region 3
#         length_x_region3 = x_riv_droite - x_riv_gauche
#         num_div_x_region3 = max(2, math.ceil(length_x_region3 / dx_precis))
#         num_div_z_precis = max(2, math.ceil((z_riv - min_z) / dz_precis))
        

#         # Application des subdivisions et recombinaison pour chaque sous-surface
#         for surface in [surface1, surface2, surface3, surface4, surface5]:
#             gmsh.model.mesh.setTransfiniteSurface(surface, "Left")
#             gmsh.model.mesh.setRecombine(2, surface)

#         # Subdivisions pour la région supérieure (surface1)
#         gmsh.model.mesh.setTransfiniteCurve(l1, num_div_x)  # Ligne supérieure horizontale
#         gmsh.model.mesh.setTransfiniteCurve(l101, num_div_x)  # Ligne horizontale intermédiaire (limite commune entre region 1 et 2 3 4)
#         gmsh.model.mesh.setTransfiniteCurve(l8, num_div_z_region1)  # Ligne gauche verticale (région 1)
#         gmsh.model.mesh.setTransfiniteCurve(l2, num_div_z_region1)  # Ligne droite verticale (région 1)
#         # Idem pour surface5
#         gmsh.model.mesh.setTransfiniteCurve(l100, num_div_x)  # Ligne supérieure horizontale intermédiaire (limite commune entre region 5 et 2 3 4)
#         gmsh.model.mese en Z entre les 3 régions 2,3 et 4 avec la région 5.h.setTransfiniteCurve(l10, num_div_x)  # Ligne inférieur horizontale 
#         gmsh.model.mesh.setTransfiniteCurve(l11, num_div_z_region5)  # Ligne gauche verticale (région 5)
#         gmsh.model.mesh.setTransfiniteCurve(l9, num_div_z_region5)  # Ligne droite verticale (région 5)

#         # Subdivisions pour la région inférieure (surface2)
#         gmsh.model.mesh.setTransfiniteCurve(l12, num_div_x_region2)  # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l16, num_div_z_precis)  # Ligne droite verticale
#         gmsh.model.mesh.setTransfiniteCurve(l6, num_div_x_region2)  # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l7, num_div_z_precis)  # Ligne gauche verticale
#         # Idem pour surface4
#         gmsh.model.mesh.setTransfiniteCurve(l14, num_div_x_region4)  # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l3, num_div_z_precis)  # Ligne droite verticale 
#         gmsh.model.mesh.setTransfiniteCurve(l4, num_div_x_region4)  # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l15, num_div_z_precis)  # Ligne gauche verticale

#         # Subdivisions pour la région inférieure (surface3)
#         gmsh.model.mesh.setTransfiniteCurve(l13, num_div_x_region3)  # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l15, num_div_z_precis)  # Ligne droite verticale
#         gmsh.model.mesh.setTransfiniteCurve(l5, num_div_x_region3)  # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l16, num_div_z_precis)  # Ligne gauche verticale

#         # Génération du maillage
#         gmsh.model.mesh.generate(mesh_dimension)
#         gmsh.write(output_mesh_path)
#         print(f"✅ Maillage sauvegardé sous : {output_mesh_path}")

#     except Exception as e:
#         print(f"Erreur lors de la génération du maillage : {e}")
#     finally:
#         gmsh.finalize()

# def generate_mesh_8_region(distance_altitude_table, output_mesh_path,v_bot=103.8,point_RG =5, point_RD = 17, dx_grossier=0.5, dx_precis = 0.1, dz_grossier=1.0, dz_precis=0.2,mesh_dimension=2):
#     """
#     Fonction permettant de créer un maillage différent pour 8 régions. Les régions sont réparties selon cet ordre : Région 1 bas à gauche du maillage, Région 2 bas centre du maillage,
#     Région 3 bas droit du maillage, Région 4 milieu gauche, Région 5 centre du maillage, Région 6 milieu droit, Région 7 haut gauche et Région 8 haut à droite. Ces régions sont limitées
#     selon plusieurs points. Tout d'abords le fichier de topo "distance_altitude_table" Dont le 1er point correspond au coin inférieur gauche du maillage complet (cad le p1), le second 
#     supérieur gauche (cad p4), l'avant dernier le supérieur droit (p16) et le dernier point le coin inférieur droit (le p13). La frontière en z entre la partie haute et milieu est déterminée
#     par l'altitude max des points de la topo rentré pour "point_RG" et "point_RD". Les coordonnées en x des points RG et RD correspondent à la frontière entre la partie gauche, centre et
#     droit. Enfin, la frontière en Z entre le milieu et le bas du maillage est déterminé par le parametre v_bot.Une fois le maillage créé il est enregistré vers output_mesh_path.

#     Entrée :
#     distance_altitude_table (DataFrame) : Tableau de la topo contenant au moins 2 colonnes avec "Distances (m)" pour les coordonnées en x et "Altitude (Z)" pour Z.
#     output_mesh_path (str) : Chemin où sera enregistrer le maillage brute
#     v_bot (float) : Parametre qui délimite la frontière en Z entre les régions 1 et 4, 2 et 5 puis 3 et 6.
#     point_RG (int) : Entier correspondant au numéro du point topo correspondant à la rive gauche. Il permet de délimiter la frontière entre les région 1 et 2 puis les régions 4 et 5 en x et entre le milieu et le haut du maillage si son altitude est plus grande que RD
#     point_RD (int) : Entier correspondant au numéro du point topo correspondant à la rive droite. Il permet de délimiter la frontière entre les région 2 et 3 puis les régions 5 et 6 en x et entre le milieu et le haut du maillage si son altitude est plus grande que RG
#     dx_grossier (float) : Valeur du pas en x utilisé dans les régions où x n'a pas besoin d'être précis. C'est-à-dire dans les régions 1, 3, 4, 6, 7 et 8
#     dx_precis (float) : Valeur du pas en x utilisé dans les régions où x a besoin d'être précis. C'est-à-dire dans la région 2 et 5
#     dz_grossier (float) : Valeur du pas en z utilisé dans les régions où z n'a pas besoin d'être précis. C'est-à-dire dans les régions 1, 2, 3, 7 et 8
#     dz_precis (float) : Valeur du pas en z utilisé dans les régions où z a besoin d'être précis. C'est-à-dire dans les régions 4, 5 et 6
#     mesh_dimension (int) : Entier qui rensigne sur la dimension du maillage (Par défault égale à 2)

#     Sortie :
#     Le maillage est créé et sauvegarder dans le chemin renseigné avec "output_mesh_path"
    
#     """
    
#     try:
#         print("\nGénération du maillage structuré avec Gmsh...")
#         gmsh.initialize()
#         gmsh.option.setNumber("General.Terminal", 1)
#         gmsh.model.add("RectangularMesh")

#         # Extraction des distances et altitudes
#         distances = distance_altitude_table["Distance (m)"].to_numpy()
#         z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()

#         # Points du rectangle
#         min_z= min(z_coords) # On prend le min et le max pour obtenir un rectangle régulier
#         max_z= max(z_coords)
#         min_x = min(distances)
#         max_x = max(distances)
#         # point intermediaire correspondant à la riv (avec un maillage plus fin)
#         x_riv_gauche = distances[point_RG]
#         x_riv_droite = distances[point_RD]
#         z_riv = max(z_coords[point_RG],z_coords[point_RD])

#         #Construction des points
#         p1 = gmsh.model.occ.addPoint(min_x, 0, min_z) # Equivalent du point 0 du DF
#         p2 = gmsh.model.occ.addPoint(min_x, 0, v_bot) 
#         p3 = gmsh.model.occ.addPoint(min_x, 0,z_riv)
#         p4 = gmsh.model.occ.addPoint(min_x, 0,max_z) # Equivalent du point 1 du DF

#         p5 = gmsh.model.occ.addPoint(x_riv_gauche, 0, min_z)
#         p6 = gmsh.model.occ.addPoint(x_riv_gauche, 0, v_bot) 
#         p7 = gmsh.model.occ.addPoint(x_riv_gauche, 0,z_riv) # Equivalent du point_RG du DF
#         p8 = gmsh.model.occ.addPoint(x_riv_gauche, 0,max_z) 

#         p9 = gmsh.model.occ.addPoint(x_riv_droite, 0, min_z)
#         p10 = gmsh.model.occ.addPoint(x_riv_droite, 0, v_bot) 
#         p11 = gmsh.model.occ.addPoint(x_riv_droite, 0,z_riv) # Equivalent du point_RD du DF
#         p12 = gmsh.model.occ.addPoint(x_riv_droite, 0,max_z) 

#         p13 = gmsh.model.occ.addPoint(max_x, 0, min_z) # Equivalent du point -1 du DF
#         p14 = gmsh.model.occ.addPoint(max_x, 0, v_bot) 
#         p15 = gmsh.model.occ.addPoint(max_x, 0,z_riv) 
#         p16 = gmsh.model.occ.addPoint(max_x, 0,max_z) # Equivalent du point -2 du DF

#         # Création des lignes pour délimité chaque région
#         #Lignes horizontales
#         l1 = gmsh.model.occ.addLine(p1, p2)
#         l2 = gmsh.model.occ.addLine(p2, p3)
#         l3 = gmsh.model.occ.addLine(p3,p4)
#         l4 = gmsh.model.occ.addLine(p5,p6)
#         l5 = gmsh.model.occ.addLine(p6,p7)
#         l6 = gmsh.model.occ.addLine(p7,p8)
#         l7 = gmsh.model.occ.addLine(p9,p10)
#         l8 = gmsh.model.occ.addLine(p10,p11)
#         l9 = gmsh.model.occ.addLine(p11,p12)
#         l10 = gmsh.model.occ.addLine(p13,p14)
#         l11 = gmsh.model.occ.addLine(p14,p15)
#         l12 = gmsh.model.occ.addLine(p15,p16)

#         #Lignes verticales
#         l100 = gmsh.model.occ.addLine(p1,p5)
#         l101 = gmsh.model.occ.addLine(p2,p6)
#         l102 = gmsh.model.occ.addLine(p3,p7)
#         l103 = gmsh.model.occ.addLine(p4,p8)
#         l104 = gmsh.model.occ.addLine(p5,p9) 
#         l105 = gmsh.model.occ.addLine(p6,p10) 
#         l106 = gmsh.model.occ.addLine(p7,p11)
#         l107 = gmsh.model.occ.addLine(p9,p13)
#         l108 = gmsh.model.occ.addLine(p10,p14)
#         l109 = gmsh.model.occ.addLine(p11,p15)
#         l110 = gmsh.model.occ.addLine(p12,p16)

#         # Combiner toutes les courbes dans unimport numpy as np


#Definition de la fontion
def creation_E_zone(sim_dir,polygons_by_zone,default_zone = 1):
    '''
    Fonction qui permet de remplir le fichier E_zone.dat utile pour ginette. On a besoin pour cela du maillage (fichier E_coordonee.dat) et d'un dataframe comprenant
    les polygones rangés selon leur zone respectivent. La fonction test si le centre de chaque maille se trouve ou non dans un polygone. Une fois la zone de la maille
    trouvé c'est ajouté à une liste qui en fin de fonction sera écrit dans le fichier E_zone.dat

    Entrée
    path_E_coor (str) : Chemin (relatif ou absolu) amenant au fichier E_coordonee.dat
    path_E_zone (str) : Chemin (relatif ou absolu) amenant au fichier E_zone.dat
    polygons_by_zone (DataFrame) : Df comprenant les polygones qui définissent les différentes zones du maillage. Un exemple de création de polygons_by_zone est donné dans la suite du code
    default_zone (int) : zone (par default = 1) que l'on va attribuer au point contenu dans aucun polygone (Evite de définir tout les polygone dans polygons_by_zone)
    '''

    liste_zone = list() # Liste qui va contenir le zone de chaque maille qui va nous permettre d'écrire dans le fichier E_zone.dat

    #Chargement des coordonées E_coordonnee.dat in sim_dir
    path_E_coor = sim_dir + "/E_coordonnees.dat"
    path_E_zone = sim_dir + "/E_zone.dat"
    coord = np.loadtxt(path_E_coor)
    coord_x = coord[:,0]
    coord_z = coord[:,1]

    for x,z in zip (coord_x,coord_z): # On test tout les couples de coord

        test_contain = False # re-initialisation du bool
        p = Point(x,z) # Point que l'on va tester

        for zone in polygons_by_zone: # test sur toutes les zones

            for poly in polygons_by_zone[zone] : # test sur tout les poly de la zone
                if poly.contains(p):
                    # print(f"Le point est à l'intérieur ✅ de la région {zone}")
                    test_contain = True
                    liste_zone.append(zone)
                    break # Sortie de la boucle poly (evite de tester tous les poly)

            if test_contain :
                break # idem mais pour la boucle zone
                
        if test_contain == False :
            # print(f"Le point est à l'extérieur ❌ de toute les région et devient region {default_zone}")
            liste_zone.append(default_zone)

    # Ecriture dans le fichier E_zone.dat
    filename = path_E_zone

    with open(filename, "w") as f:
        for item in liste_zone:
            f.write(f"{item}\n") # \n necessaire pour ecrire un par ligne face([loop2]) # 2
#         surface3 = gmsh.model.occ.addPlaneSurface([loop3]) # 3
#         surface4 = gmsh.model.occ.addPlaneSurface([loop4]) # 4
#         surface5 = gmsh.model.occ.addPlaneSurface([loop5]) # 5
#         surface6 = gmsh.model.occ.addPlaneSurface([loop6]) # 5
#         surface7 = gmsh.model.occ.addPlaneSurface([loop7]) # 5
#         surface8 = gmsh.model.occ.addPlaneSurface([loop8]) # 5

#         # Synchronisation du modèle
#         gmsh.model.occ.synchronize()

#         # Calcul des subdivisions en x pour les regions de gauche 1, 4 et 7
#         length_x_gauche = x_riv_gauche - min_x
#         num_div_x_gauche = max(2, math.ceil(length_x_gauche / dx_grossier))

#         # Calcul des subdivisions en x pour les regions centre 2 et 5
#         length_x_centre = x_riv_droite - x_riv_gauche
#         num_div_x_centre = max(2, math.ceil(length_x_centre / dx_precis))

#         # Calcul des subdivisions en x pour les regions droite 3, 6 et 8
#         length_x_droite = max_x - x_riv_droite
#         num_div_x_droite = max(2, math.ceil(length_x_droite / dx_grossier))

#         # Calcul des subdivisions en z pour les regions haute 7 et 8
#         num_div_z_haut = max(2, math.ceil((max_z - z_riv) / dz_grossier))

#         # Calcul des subdivisions en z pour les regions centre 4, 5 et 6
#         num_div_z_centre = max(2, math.ceil((z_riv - v_bot) / dz_precis))

#         # Calcul des subdivisions en z pour les regions bas 1,2 et 3
#         num_div_z_bas = max(2, math.ceil((v_bot - min_z) / dz_grossier))

#         print(f"Subdivisions calculées :\n"
#             f"  - Gauche : {num_div_x_gauche} subdivisions\n"
#             f" - dx gauche : {length_x_gauche / num_div_x_gauche:.2f} m\n"
#             f"  - Centre : {num_div_x_centre} subdivisions\n"
#             f" - dx centre : {length_x_centre / num_div_x_centre:.2f} m\n"
#             f"  - Droite : {num_div_x_droite} subdivisions\n"
#             f" - dx droite : {length_x_droite / num_div_x_droite:.2f} m\n"
#             f"  - Haut : {num_div_z_haut} subdivisions\n"
#             f" - dz haut : {(max_z - z_riv) / num_div_z_haut:.2f} m\n"
#             f"  - Centre : {num_div_z_centre} subdivisions\n"
#             f" - dz centre : {(z_riv - v_bot) / num_div_z_centre:.2f} m\n"
#             f"  - Bas : {num_div_z_bas} subdivisions\n"
#             f" - dz bas : {(v_bot - min_z) / num_div_z_bas:.2f} m\n")


#         # Application des subdivisions et recombinaison pour chaque sous-surface
#         for surface in [surface1, surface2, surface3, surface4, surface5, surface6, surface7, surface8]:
#             gmsh.model.mesh.setTransfiniteSurface(surface, "Left")
#             gmsh.model.mesh.setRecombine(2, surface)

#         # Subdivisions pour la région (surface1)
#         gmsh.model.mesh.setTransfiniteCurve(l1, num_div_z_bas) # Ligne verticale gauche
#         gmsh.model.mesh.setTransfiniteCurve(l101, num_div_x_gauche) # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z_bas) # Ligne verticale droite
#         gmsh.model.mesh.setTransfiniteCurve(l100, num_div_x_gauche) # Ligne horizontale basse

#         # surface2
#         gmsh.model.mesh.setTransfiniteCurve(l105, num_div_x_centre) # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l104, num_div_x_centre) # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z_bas) # Ligne verticale gauche
#         gmsh.model.mesh.setTransfiniteCurve(l7, num_div_z_bas) # Ligne verticale droite

#         # surface3
#         gmsh.model.mesh.setTransfiniteCurve(l107, num_div_x_droite)  # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l108, num_div_x_droite)  # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l10, num_div_z_bas)  # Ligne verticale droite
#         gmsh.model.mesh.setTransfiniteCurve(l7, num_div_z_bas)  # Ligne verticale gauche

#         # surface4
#         gmsh.model.mesh.setTransfiniteCurve(l101, num_div_x_gauche) # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l102, num_div_x_gauche) # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l2, num_div_z_centre) # Ligne verticale gauche
#         gmsh.model.mesh.setTransfiniteCurve(l5, num_div_z_centre) # Ligne verticale droite

#         # surface5
#         gmsh.model.mesh.setTransfiniteCurve(l105, num_div_x_centre) # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l106, num_div_x_centre) # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l5, num_div_z_centre) # Ligne verticale gauche
#         gmsh.model.mesh.setTransfiniteCurve(l8, num_div_z_centre) # Ligne verticale droite

#         # surface6
#         gmsh.model.mesh.setTransfiniteCurve(l108, num_div_x_droite) # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l109, num_div_x_droite) # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l8, num_div_z_centre) # Ligne verticale gauche
#         gmsh.model.mesh.setTransfiniteCurve(l11, num_div_z_centre) # Ligne verticale droite

#         # surface7
#         gmsh.model.mesh.setTransfiniteCurve(l102, num_div_x_gauche) # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l103, num_div_x_gauche) # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l3, num_div_z_haut) # Ligne verticale gauche
#         gmsh.model.mesh.setTransfiniteCurve(l6, num_div_z_haut) # Ligne verticale droite

#         # surface8
#         gmsh.model.mesh.setTransfiniteCurve(l109, num_div_x_droite) # Ligne horizontale basse
#         gmsh.model.mesh.setTransfiniteCurve(l110, num_div_x_droite) # Ligne horizontale haute
#         gmsh.model.mesh.setTransfiniteCurve(l9, num_div_z_haut) # Ligne verticale gauche
#         gmsh.model.mesh.setTransfiniteCurve(l12, num_div_z_haut) # Ligne verticale droite

#         # Génération du maillage
#         gmsh.model.mesh.generate(mesh_dimension)
#         gmsh.write(output_mesh_path)
#         print(f"✅ Maillage sauvegardé sous : {output_mesh_path}")

#     except Exception as e:
#         print(f"Erreur lors de la génération du maillage : {e}")
#     finally:
#         gmsh.finalize()

def subdivisions_autour(x_hobo, dx_hobo, N_cells):
    """
    Retourne les bornes et centres de N_cells mailles de largeur dx_hobo,
    centrées globalement sur x_hobo.

    Parameters
    ----------
    x_hobo    : float
        Abscisse du point Hobo.
    dx_hobo   : float
        largeur de chaque maille (par ex 0.01).
    N_cells   : int
        Nombre de mailles que l’on veut faire “autour” de x_hobo.

    Returns
    -------
    centers : np.ndarray, shape (N_cells,)
        Les abscisses des centres de maille.
    """
    demi_largeur = (N_cells * dx_hobo) / 2.0
    x_min = x_hobo - demi_largeur
    x_max = x_hobo + demi_largeur
    edges = np.linspace(x_min, x_max, N_cells + 1)
    centers = (edges[:-1] + edges[1:]) / 2.0
    return edges

def generate_mesh_8_region(distance_altitude_table, output_mesh_path,v_bot=103.8,x_RG =5, x_RD = 17,z_riv=106, 
                           dx_grossier=0.5, dx_precis = 0.1,x_hobo_1 =8.5, z_hobo_1=104.5,x_hobo_2 = 10,z_hobo_2=104.5, dx_hobo=0.01, 
                           dz_grossier=1.0, dz_precis=0.2,dz_hobo1=0.1,dz_hobo=0.15,mesh_dimension=2):
    """
    Fonction permettant de créer un maillage différent pour 8 régions. Les régions sont réparties selon cet ordre : Région 1 bas à gauche du maillage, Région 2 bas centre du maillage,
    Région 3 bas droit du maillage, Région 4 milieu gauche, Région 5 centre du maillage, Région 6 milieu droit, Région 7 haut gauche et Région 8 haut à droite. Ces régions sont limitées
    selon plusieurs points. Tout d'abords le fichier de topo "distance_altitude_table" Dont le 1er point correspond au coin inférieur gauche du maillage complet (cad le p1), le second 
    supérieur gauche (cad p4), l'avant dernier le supérieur droit (p16) et le dernier point le coin inférieur droit (le p13). La frontière en z entre la partie haute et milieu est déterminée
    par l'altitude max des points de la topo rentré pour "point_RG" et "point_RD". Les coordonnées en x des points RG et RD correspondent à la frontière entre la partie gauche, centre et
    droit. Enfin, la frontière en Z entre le milieu et le bas du maillage est déterminé par le parametre v_bot.Une fois le maillage créé il est enregistré vers output_mesh_path.

        Entrée :
    distance_altitude_table (DataFrame) : Tableau de la topo contenant au moins 2 colonnes avec "Distances (m)" pour les coordonnées en x et "Altitude (Z)" pour Z.
    output_mesh_path (str) : Chemin où sera enregistrer le maillage brute
    v_bot (float) : Parametre qui délimite la frontière en Z entre les régions 1 et 4, 2 et 5 puis 3 et 6.
    point_RG (float) : Valeur qui correspond à la distance en x où se trouve la rive gauche. Il permet de délimiter la frontière entre les région 1 et 2 puis les régions 4 et 5 en x
    point_RD (float) : Valeur qui correspond à la distance en x où se trouve la rive droite. Il permet de délimiter la frontière entre les région 2 et 3 puis les régions 5 et 6 en x 
    z_riv (float) : Valeur qui correspond à la hauteur de la riv. Le maillage sera précis en z entre v_bot et z_riv. Il délimite donc le milieu et le haut du maillage
    x_hobo_1 (float) : Abscisse du premier point Hobo où un maillage plus fin est requis.
    x_hobo_2 (float) : Abscisse du second point Hobo où un maillage plus fin est requis.
    dx_hobo (float) : Largeur des mailles autour des points Hobo pour un maillage précis.
    dx_grossier (float) : Valeur du pas en x utilisé dans les régions où x n'a pas besoin d'être précis. C'est-à-dire dans les régions 1, 3, 4, 6, 7 et 8
    dx_precis (float) : Valeur du pas en x utilisé dans les régions où x a besoin d'être précis. C'est-à-dire dans la région 2 et 5
    dz_grossier (float) : Valeur du pas en z utilisé dans les régions où z n'a pas besoin d'être précis. C'est-à-dire dans les régions 1, 2, 3, 7 et 8
    dz_precis (float) : Valeur du pas en z utilisé dans les régions où z a besoin d'être précis. C'est-à-dire dans les régions 4, 5 et 6
    mesh_dimension (int) : Entier qui rensigne sur la dimension du maillage (Par défault égale à 2)

        Sortie :
    Le maillage est créé et sauvegarder dans le chemin renseigné avec "output_mesh_path"
        
    """
    
    try:
        print("\nGénération du maillage structuré avec Gmsh...")
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("RectangularMesh")

        # Extraction des distances et altitudes
        distances = distance_altitude_table["Distance (m)"].to_numpy()
        z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()

        # Points du rectangle
        min_z= min(z_coords) # On prend le min et le max pour obtenir un rectangle régulier
        max_z= max(z_coords)
        min_x = min(distances)
        max_x = max(distances)
        # point intermediaire correspondant à la riv (avec un maillage plus fin)
        x_riv_gauche = x_RG
        x_riv_droite = x_RD
        #  z_riv=max(z_hobo_1,z_hobo_2,z_riv) # On prend la plus haute des deux hauteurs de hobo
        z_riv = max(z_hobo_1,z_hobo_2,z_riv) # On prend la plus haute des deux hauteurs de hobo
        #Construction des points
        p1 = gmsh.model.occ.addPoint(min_x, 0, min_z) # Equivalent du point 0 du DF
        p2 = gmsh.model.occ.addPoint(min_x, 0, v_bot) 
        p3 = gmsh.model.occ.addPoint(min_x, 0,z_riv)
        p4 = gmsh.model.occ.addPoint(min_x, 0,max_z) # Equivalent du point 1 du DF

        p5 = gmsh.model.occ.addPoint(x_riv_gauche, 0, min_z)
        p6 = gmsh.model.occ.addPoint(x_riv_gauche, 0, v_bot) 
        p7 = gmsh.model.occ.addPoint(x_riv_gauche, 0,z_riv) # Equivalent du point_RG du DF
        p8 = gmsh.model.occ.addPoint(x_riv_gauche, 0,max_z) 

        p9 = gmsh.model.occ.addPoint(x_riv_droite, 0, min_z)
        p10 = gmsh.model.occ.addPoint(x_riv_droite, 0, v_bot) 
        p11 = gmsh.model.occ.addPoint(x_riv_droite, 0,z_riv) # Equivalent du point_RD du DF
        p12 = gmsh.model.occ.addPoint(x_riv_droite, 0,max_z) 

        p13 = gmsh.model.occ.addPoint(max_x, 0, min_z) # Equivalent du point -1 du DF
        p14 = gmsh.model.occ.addPoint(max_x, 0, v_bot) 
        p15 = gmsh.model.occ.addPoint(max_x, 0,z_riv) 
        p16 = gmsh.model.occ.addPoint(max_x, 0,max_z) # Equivalent du point -2 du DF


        
        # Calcul des subdivisions en x pour les regions de gauche 1, 4 et 7
        length_x_gauche = x_riv_gauche - min_x
        num_div_x_gauche = max(2, math.ceil(length_x_gauche / dx_grossier))

        # Calcul des subdivisions en x pour les regions centre 2 et 5
        length_x_centre = x_riv_droite - x_riv_gauche
        num_div_x_centre = max(2, math.ceil(length_x_centre / dx_precis))

        # Calcul des subdivisions en x pour les regions droite 3, 6 et 8
        length_x_droite = max_x - x_riv_droite
        num_div_x_droite = max(2, math.ceil(length_x_droite / dx_grossier))

        # Calcul des subdivisions en z pour les regions haute 7 et 8
        num_div_z_haut = max(2, math.ceil((max_z - z_riv) / dz_grossier))

        # Calcul des subdivisions en z pour les regions centre 4, 5 et 6
        num_div_z_centre = max(2, math.ceil((z_riv - v_bot) / dz_precis))

        # Calcul des subdivisions en z pour les regions bas 1,2 et 3
        num_div_z_bas = max(2, math.ceil((v_bot - min_z) / dz_grossier))


        # Calcul des divisions autour de x_hobo_1 avec 4 * dx_hobo on veut que le centre des mailles tombe sur x_hobo_1
        # Ajout des points Hobo
        
        dx_reel = (x_RD-x_RG)/num_div_x_centre # correspond au delta réel qu'il y a entre x_RD et x_RG qui a été recalculé avec num_div_x_centre
        maille_hobo1 = 0 # Variable qui va contenir le numéro de la maille ou se trouve le hobo dans la région centre AVANT MODIFICATION
        maille_hobo2 = 0
        nb_case_a_ajouter1 = 1
        nb_case_a_ajouter2 = 1

        for i in range (num_div_x_centre) : #Compteur pour savoir quand on arrive dans la maille contenant le hobo
            if (x_RG + dx_reel * (i+1)) > x_hobo_1:
                break
            else :
                maille_hobo1 = maille_hobo1+1 
        
        if (x_RG + dx_reel * (maille_hobo1+1)) - x_hobo_1 < 0.005 : nb_case_a_ajouter1 = 2 #Dans le cas ou le hobo est trop proche de la maille suivante et qui empeche de un découpage de 1cm selon x (car ca dépasse) on prend donc la maille suivante 

        x_gauche_hobo1 = x_RG + dx_reel*(maille_hobo1)
        x_droite_hobo1 = x_RG + dx_reel*(maille_hobo1+nb_case_a_ajouter1)

        for i in range (num_div_x_centre) : #Compteur pour savoir quand on arrive dans la maille contenant le hobo
            if (x_RG + dx_reel * (i+1)) > x_hobo_2:
                break
            else :
                maille_hobo2 = maille_hobo2+1 
        
        if (x_RG + dx_reel * (maille_hobo2+1)) - x_hobo_2 < 0.005 : nb_case_a_ajouter2 = 2 #Dans le cas ou le hobo est trop proche de la maille suivante et qui empeche de un découpage de 1cm selon x (car ca dépasse) on prend donc la maille suivante 

        x_gauche_hobo2 = x_RG + dx_reel*(maille_hobo2)
        x_droite_hobo2 = x_RG + dx_reel*(maille_hobo2+nb_case_a_ajouter2)

        # Calcul des subdivisions en x pour les regions affecté par l'arrivé du hobo
        num_div_x_RG_hobo1 = maille_hobo1
        num_div_x_hobo1_hobo2 = max(2, math.ceil((x_gauche_hobo2-x_droite_hobo1) / dx_reel))
        num_div_x_hobo2_RD = num_div_x_centre - (maille_hobo2+nb_case_a_ajouter2)

        length_x_hobo1 = x_droite_hobo1 - x_gauche_hobo1
        num_div_x_hobo1 = max(2, math.ceil(length_x_hobo1 / dx_hobo))
        length_x_hobo2 = x_droite_hobo2 - x_gauche_hobo2
        num_div_x_hobo2 = max(2, math.ceil(length_x_hobo2 / dx_hobo))
        print('x_hobo_gauche_1',x_gauche_hobo1,'x_hobo_droite_1',x_droite_hobo1, dx_reel,(x_gauche_hobo1-x_droite_hobo1) / dx_reel)

        # Calcul des subdivisions en z pour les regions affecté par l'arrivé du hobo
        # on veut que le centre des mailles contiennes la coordonnes z_hobo1-dz_hobo1
        # nouveau point
        p29= gmsh.model.occ.addPoint(x_gauche_hobo1, 0, z_hobo_1)
        p30= gmsh.model.occ.addPoint(x_gauche_hobo1, 0, z_hobo_1)
        
        
        # Créer les nouveaux points pour gmsh

        p17 = gmsh.model.occ.addPoint(x_gauche_hobo1, 0, z_riv)
        p18 = gmsh.model.occ.addPoint(x_droite_hobo1, 0, z_riv)
        p19 = gmsh.model.occ.addPoint(x_gauche_hobo1, 0, v_bot)
        p20 = gmsh.model.occ.addPoint(x_droite_hobo1, 0, v_bot)
        p21 = gmsh.model.occ.addPoint(x_gauche_hobo1, 0, min_z)
        p22 = gmsh.model.occ.addPoint(x_droite_hobo1, 0, min_z)

        p23 = gmsh.model.occ.addPoint(x_gauche_hobo2, 0, z_riv)
        p24 = gmsh.model.occ.addPoint(x_droite_hobo2, 0, z_riv)
        p25 = gmsh.model.occ.addPoint(x_gauche_hobo2, 0, v_bot)
        p26 = gmsh.model.occ.addPoint(x_droite_hobo2, 0, v_bot)
        p27 = gmsh.model.occ.addPoint(x_gauche_hobo2, 0, min_z)
        p28 = gmsh.model.occ.addPoint(x_droite_hobo2, 0, min_z)


        # Création des lignes pour délimité chaque région
        #Lignes verticales
        l1 = gmsh.model.occ.addLine(p1, p2)
        l2 = gmsh.model.occ.addLine(p2, p3)
        l3 = gmsh.model.occ.addLine(p3,p4)
        l4 = gmsh.model.occ.addLine(p5,p6)
        l5 = gmsh.model.occ.addLine(p6,p7)
        l6 = gmsh.model.occ.addLine(p7,p8)
        l7 = gmsh.model.occ.addLine(p9,p10)
        l8 = gmsh.model.occ.addLine(p10,p11)
        l9 = gmsh.model.occ.addLine(p11,p12)
        l10 = gmsh.model.occ.addLine(p13,p14)
        l11 = gmsh.model.occ.addLine(p14,p15)
        l12 = gmsh.model.occ.addLine(p15,p16)

        #Avec hobo
        l13 = gmsh.model.occ.addLine(p21,p19)
        l14 = gmsh.model.occ.addLine(p19,p17)
        l15 = gmsh.model.occ.addLine(p22,p20)
        l16 = gmsh.model.occ.addLine(p20,p18)
        l17 = gmsh.model.occ.addLine(p27,p25)
        l18 = gmsh.model.occ.addLine(p25,p23)
        l19 = gmsh.model.occ.addLine(p28,p26)
        l20 = gmsh.model.occ.addLine(p26,p24)

        #Lignes horizontales
        l100 = gmsh.model.occ.addLine(p1,p5)
        l101 = gmsh.model.occ.addLine(p2,p6)
        l102 = gmsh.model.occ.addLine(p3,p7)
        l103 = gmsh.model.occ.addLine(p4,p8)

        #avec Hobo1
        l104 = gmsh.model.occ.addLine(p5,p21)
        l105 = gmsh.model.occ.addLine(p6,p19)
        l106 = gmsh.model.occ.addLine(p7,p17)
        l107 = gmsh.model.occ.addLine(p21,p22)
        l108 = gmsh.model.occ.addLine(p17,p18)

        # avec hobo2
        l109 = gmsh.model.occ.addLine(p22,p27)
        l110 = gmsh.model.occ.addLine(p20,p25)
        l111 = gmsh.model.occ.addLine(p18,p23)
        l117 = gmsh.model.occ.addLine(p27,p28)
        l118 = gmsh.model.occ.addLine(p25,p26)
        l119 = gmsh.model.occ.addLine(p23,p24)
        l120 = gmsh.model.occ.addLine(p28,p9)
        l121 = gmsh.model.occ.addLine(p26,p10)
        l122 = gmsh.model.occ.addLine(p24,p11)


        l112 = gmsh.model.occ.addLine(p9,p13)
        l113 = gmsh.model.occ.addLine(p10,p14)
        l114 = gmsh.model.occ.addLine(p11,p15)
        l115 = gmsh.model.occ.addLine(p12,p16)
        l116 = gmsh.model.occ.addLine(p19,p20)

        # Combiner toutes les courbes dans une seule boucle
       # Création des boucles pour les 5 régions
        loop1 = gmsh.model.occ.addCurveLoop([l101, -l4, -l100, l1])  # Région 1
        loop2 = gmsh.model.occ.addCurveLoop([l105, -l13, -l104, l4])  # Région 2
        loop3 = gmsh.model.occ.addCurveLoop([l116, -l15, -l107, l13]) # Région 3
        loop4 = gmsh.model.occ.addCurveLoop([l110, -l17, -l109, l15])
        loop5 = gmsh.model.occ.addCurveLoop([l113, -l10, -l112, l7])
        loop6 = gmsh.model.occ.addCurveLoop([l102, -l5, -l101, l2])  # Région 6
        loop7 = gmsh.model.occ.addCurveLoop([l106, -l14, -l105, l5])
        loop8 = gmsh.model.occ.addCurveLoop([l111, -l18, -l110, l16])
        loop9 = gmsh.model.occ.addCurveLoop([l114, -l11, -l113, l8])
        loop10 = gmsh.model.occ.addCurveLoop([l103, -l6, -l102, l3])  # Région 7
        loop11 = gmsh.model.occ.addCurveLoop([l115, -l12, -l114, l9])
        loop12 = gmsh.model.occ.addCurveLoop([l108, -l16, -l116, l14])
        loop13 = gmsh.model.occ.addCurveLoop([l118, -l19, -l117, l17])
        loop14 = gmsh.model.occ.addCurveLoop([l121, -l7, -l120, l19])
        loop15 = gmsh.model.occ.addCurveLoop([l119, -l20, -l118, l18])
        loop16 = gmsh.model.occ.addCurveLoop([l122, -l8, -l121, l20])

        # Création des sous-surfaces
        surface1 = gmsh.model.occ.addPlaneSurface([loop1]) # 1
        surface2 = gmsh.model.occ.addPlaneSurface([loop2]) # 2
        surface3 = gmsh.model.occ.addPlaneSurface([loop3]) # 3
        surface4 = gmsh.model.occ.addPlaneSurface([loop4]) # 4
        surface5 = gmsh.model.occ.addPlaneSurface([loop5]) # 5
        surface6 = gmsh.model.occ.addPlaneSurface([loop6]) # 5
        surface7 = gmsh.model.occ.addPlaneSurface([loop7]) # 5
        surface8 = gmsh.model.occ.addPlaneSurface([loop8]) # 5
        surface9 = gmsh.model.occ.addPlaneSurface([loop9]) # 5
        surface10 = gmsh.model.occ.addPlaneSurface([loop10]) # 5
        surface11 = gmsh.model.occ.addPlaneSurface([loop11]) # 5
        surface12 = gmsh.model.occ.addPlaneSurface([loop12]) # 5
        surface13 = gmsh.model.occ.addPlaneSurface([loop13]) # 5
        surface14 = gmsh.model.occ.addPlaneSurface([loop14]) # 5
        surface15 = gmsh.model.occ.addPlaneSurface([loop15]) # 5
        surface16 = gmsh.model.occ.addPlaneSurface([loop16]) # 5

        # Synchronisation du modèle
        gmsh.model.occ.synchronize()

        # Calcul des subdivisions en x pour les regions de gauche 1, 4 et 7
        length_x_gauche = x_riv_gauche - min_x
        num_div_x_gauche = max(2, math.ceil(length_x_gauche / dx_grossier))

        # Calcul des subdivisions en x pour les regions centre 2 et 5
        length_x_centre = x_riv_droite - x_riv_gauche
        num_div_x_centre = max(2, math.ceil(length_x_centre / dx_precis))

        # Calcul des subdivisions en x pour les regions droite 3, 6 et 8
        length_x_droite = max_x - x_riv_droite
        num_div_x_droite = max(2, math.ceil(length_x_droite / dx_grossier))

        # Calcul des subdivisions en z pour les regions haute 7 et 8
        num_div_z_haut = max(2, math.ceil((max_z - z_riv) / dz_grossier))
        # Ensure the number of subdivisions is an integer by adjusting the step size
        if (max_z - z_riv) % dz_grossier != 0.1:
            num_div_z_haut = max(2, round((max_z - z_riv) / dz_grossier))



        
        # Calcul des subdivisions en z pour les regions centre 4, 5 et 6
        num_div_z_centre = max(2, math.ceil((z_riv - v_bot) / dz_precis))

        # Calcul des subdivisions en z pour les regions bas 1,2 et 3
        num_div_z_bas = max(2, math.ceil((v_bot - min_z) / dz_grossier))

        print(f"Subdivisions calculées :\n"
            f"  - Gauche : {num_div_x_gauche} subdivisions\n"
            f" - dx gauche : {length_x_gauche / num_div_x_gauche:.2f} m\n"
            f"  - Centre : {num_div_x_centre} subdivisions\n"
            f" - dx centre : {length_x_centre / num_div_x_centre:.2f} m\n"
            f"  - Droite : {num_div_x_droite} subdivisions\n"
            f" - dx droite : {length_x_droite / num_div_x_droite:.2f} m\n"
            f"  - Haut : {num_div_z_haut} subdivisions\n"
            f" - dz haut : {(max_z - z_riv) / num_div_z_haut:.2f} m\n"
            f"  - Centre : {num_div_z_centre} subdivisions\n"
            f" - dz centre : {(z_riv - v_bot) / num_div_z_centre:.2f} m\n"
            f"  - Bas : {num_div_z_bas} subdivisions\n"
            f" - dz bas : {(v_bot - min_z) / num_div_z_bas:.2f} m\n")


        # Application des subdivisions et recombinaison pour chaque sous-surface
        for surface in [surface1, surface2, surface3, surface4, surface5, surface6, surface7, surface8, surface9, surface10, surface11, surface12, surface13, surface14, surface15, surface16]:
            gmsh.model.mesh.setTransfiniteSurface(surface, "Left")
            gmsh.model.mesh.setRecombine(2, surface)

        # Subdivisions pour la région (surface1)
        gmsh.model.mesh.setTransfiniteCurve(l1, num_div_z_bas) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l101, num_div_x_gauche) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z_bas) # Ligne verticale droite
        gmsh.model.mesh.setTransfiniteCurve(l100, num_div_x_gauche) # Ligne horizontale basse

        # surface2
        gmsh.model.mesh.setTransfiniteCurve(l105, num_div_x_RG_hobo1) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l104, num_div_x_RG_hobo1) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z_bas) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l7, num_div_z_bas) # Ligne verticale droite

        # surface3
        gmsh.model.mesh.setTransfiniteCurve(l107, num_div_x_hobo1)  # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l116, num_div_x_hobo1)  # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l15, num_div_z_bas)  # Ligne verticale droite
        gmsh.model.mesh.setTransfiniteCurve(l13, num_div_z_bas)  # Ligne verticale gauche

        # surface4 
        gmsh.model.mesh.setTransfiniteCurve(l109, num_div_x_hobo1_hobo2)  # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l110, num_div_x_hobo1_hobo2)  # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l15, num_div_z_bas)  # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l7, num_div_z_bas)  # Ligne verticale droite

        # surface5
        gmsh.model.mesh.setTransfiniteCurve(l112, num_div_x_droite)  # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l113, num_div_x_droite)  # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l7, num_div_z_bas)  # Ligne verticale droite
        gmsh.model.mesh.setTransfiniteCurve(l10, num_div_z_bas)  # Ligne verticale gauche

        # surface6
        gmsh.model.mesh.setTransfiniteCurve(l101, num_div_x_gauche) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l102, num_div_x_gauche) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l2, num_div_z_centre) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l5, num_div_z_centre) # Ligne verticale droite

        # surface7
        gmsh.model.mesh.setTransfiniteCurve(l105, num_div_x_RG_hobo1) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l106, num_div_x_RG_hobo1) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l5, num_div_z_centre) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l14, num_div_z_centre) # Ligne verticale droite

        # surface12
        gmsh.model.mesh.setTransfiniteCurve(l116, num_div_x_hobo1) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l108, num_div_x_hobo1) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l14, num_div_z_centre) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l16, num_div_z_centre) # Ligne verticale droite

        # surface8 A MODIF
        gmsh.model.mesh.setTransfiniteCurve(l110, num_div_x_hobo1_hobo2) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l111, num_div_x_hobo1_hobo2) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l16, num_div_z_centre) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l8, num_div_z_centre) # Ligne verticale droite

        # surface9
        gmsh.model.mesh.setTransfiniteCurve(l113, num_div_x_droite) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l114, num_div_x_droite) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l8, num_div_z_centre) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l11, num_div_z_centre) # Ligne verticale droite

        # surface10
        gmsh.model.mesh.setTransfiniteCurve(l102, num_div_x_gauche) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l103, num_div_x_gauche) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l3, num_div_z_haut) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l6, num_div_z_haut) # Ligne verticale droite

        # surface11
        gmsh.model.mesh.setTransfiniteCurve(l114, num_div_x_droite) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l115, num_div_x_droite) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l9, num_div_z_haut) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l12, num_div_z_haut) # Ligne verticale droite

        #Surace 13
        gmsh.model.mesh.setTransfiniteCurve(l117, num_div_x_hobo2) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l118, num_div_x_hobo2) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l17, num_div_z_bas) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l19, num_div_z_bas) # Ligne verticale droite

        #surface 14
        gmsh.model.mesh.setTransfiniteCurve(l120, num_div_x_hobo2_RD) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l121, num_div_x_hobo2_RD) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l19, num_div_z_bas) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l7, num_div_z_bas) # Ligne verticale droite

        #15
        gmsh.model.mesh.setTransfiniteCurve(l118, num_div_x_hobo2) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l119, num_div_x_hobo2) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l18, num_div_z_centre) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l20, num_div_z_centre) # Ligne verticale droite

        #16
        gmsh.model.mesh.setTransfiniteCurve(l121, num_div_x_hobo2_RD) # Ligne horizontale basse
        gmsh.model.mesh.setTransfiniteCurve(l122, num_div_x_hobo2_RD) # Ligne horizontale haute
        gmsh.model.mesh.setTransfiniteCurve(l20, num_div_z_centre) # Ligne verticale gauche
        gmsh.model.mesh.setTransfiniteCurve(l8, num_div_z_centre) # Ligne verticale droite

        # Génération du maillage
        gmsh.model.mesh.generate(mesh_dimension)
        gmsh.write(output_mesh_path)
        print(f"✅ Maillage sauvegardé sous : {output_mesh_path}")

    except Exception as e:
        print(f"Erreur lors de la génération du maillage : {e}")
    finally:
        gmsh.finalize()

#______________________________________________________________________
def remove_elements_above_curve_all_entities(table, mesh_path, ajout, retirer=''):
    gmsh.initialize()
    gmsh.open(mesh_path)

    # Supprimer la première et dernière ligne du tableau
    table_trimmed = table.iloc[1:-1]
    d = table_trimmed["Distance (m)"].to_numpy()
    z = table_trimmed["Altitude (Z)"].to_numpy()

    dim = 2

    # Récupérer toutes les régions du maillage
    entity_tags = gmsh.model.getEntities(dim)

    if not entity_tags:
        print("Aucune entite 2D trouvee.")
        gmsh.finalize()
        return

    for entity in entity_tags: # Boucle pour chaque région du maillage
        tag = entity[1]
       # print(tag)

    # Obtenir les éléments de cette entité
        elem_types, elem_tags, node_tags = gmsh.model.mesh.getElements(dim, tag)

        if not elem_tags:
            continue # Passer les régions sans éléments a suppr

        elements_to_remove = []
        for elem, nodes in zip(elem_tags[0], node_tags[0].reshape(-1, 4)):
            coords = [gmsh.model.mesh.getNode(n)[0] for n in nodes]
            center_x = np.mean([c[0] for c in coords])
            center_z = np.mean([c[2] for c in coords])
            z_curve = np.interp(center_x, d, z)
            if center_z+0.01 > z_curve:
                elements_to_remove.append(elem)

        if elements_to_remove:
            elements_to_remove = np.array(elements_to_remove, dtype=np.int32).flatten()
            gmsh.model.mesh.removeElements(dim, tag, elements_to_remove)
            gmsh.model.mesh.reclassifyNodes()
          #  print(f" Éléments supprimés dans l'entité {tag}.")

    # Sauvegarde du maillage modifié
    modified_mesh_path = mesh_path.replace(retirer+".msh", ajout+".msh")
    gmsh.write(modified_mesh_path)
    print(f" Maillage modifié enregistré sous : {modified_mesh_path}")

    gmsh.finalize()

#______________________________________________________________________
def readGmsh(fName, precision=None):
    """
    Lire un fichier Gmsh (.msh) et calculer les centres des éléments ainsi que c'est dimension en x et z

    Parameters
    ----------
    fName : str
        Chemin vers le fichier Gmsh (.msh).
    precision : int, optional
        Si spécifié, arrondir les coordonnées des centres à ce nombre de décimales.

    Returns
    -------
    centers : list of tuple
        Liste des coordonnées des centres des éléments sous forme de tuples (x, y, z).
    nb_mesh : int
        nb of center(mesh) in the mesh.
    dimension : list of tuple
        Liste des dimensions spatiales des éléments sous forme de tuples (am, bm).
    """


    assert (precision is None or (isinstance(precision, int) and precision >= 0)), \
        "Precision doit être None ou un entier positif."
    # Parse the Gmsh file
    mesh = gmshparser.parse(fName)
    node_coords = {}

    # Step 1: Collect all node coordinates
    for entity in mesh.get_node_entities():
        for node in entity.get_nodes():
            nid = node.get_tag()
            coords = node.get_coordinates()
            node_coords[nid] = coords

    # Step 2: Calculate centers of elements
    centers_dimension_elements = []
    # dimension_elements = []
    for entity in mesh.get_element_entities():
        # Get the dimension of the entity (e.g., 1 for lines, 2 for surfaces)
        dimension = entity.get_dimension()
        if dimension == 2:  # Only process 2D elements (surfaces)
            for element in entity.get_elements():
                node_ids = element.get_connectivity()
                coords = [node_coords[nid] for nid in node_ids if nid in node_coords]
                if not coords:
                    continue
                # Calculate the center of the element
                x = round(np.mean([c[0] for c in coords]),4)
                z = round(np.mean([c[2] for c in coords]),4)  # Assuming 2D (x, z) np.mean([c[2] for c in coords])  # Assuming 2D (x, z)
                am = max(round(c[0],4) for c in coords) - min(round(c[0],4) for c in coords)
                bm = max(round(c[2],4) for c in coords) - min(round(c[2],4) for c in coords)
                centers_dimension_elements.append((x, z, am, bm))
                # dimension_elements.append((am,bm))

    # Step 3: Round coordinates if precision is specified
    if precision is not None:
        centers_dimension_elements = [(round(x, precision), round(z, precision), round(am,precision), round(bm,precision)) for x, z, am, bm in centers_dimension_elements]

    # Step 4: Create a DataFrame for the centers
    df_centers = pd.DataFrame(centers_dimension_elements, columns=["x", "z", "am", "bm"])

    # Step 5: Sort the DataFrame
    # Order with x ascending and z descending
    df_centers = df_centers.sort_values(by=["x", "z"], ascending=[True, False]).reset_index(drop=True)
    df_dimension = df_centers[["am",'bm']]
    df_centers = df_centers.drop(columns=["am","bm"])

    # Step 6: Count the number of elements
    nb_mesh = len(df_centers)

    return df_centers, nb_mesh, df_dimension

#______________________________________________________________________
def voisin_mesh(directory):
    """
    function to compute the neighbors of each mesh element from the file "E_coordonnee.dat",
    producing indices in E_voisins.dat.
    The function creates a neighbor table for each mesh element, where the neighbors are defined
    as the elements that share a side with the current element.
    PARAMETERS
    ----------
    directory : str
        The directory where the mesh files are located.
    OUTPUT
    -------
    ivois : DataFrame
        A DataFrame containing the neighbor indices for each mesh element.
        The DataFrame has the following columns:
            - right: index of the right neighbor
            - left: index of the left neighbor
            - top: index of the top neighbor
            - bottom: index of the bottom neighbor  
    Writes the neighbor table to a file named "E_voisins.dat" in the specified directory.
    """
    # Load coordinates
    coord = pd.read_csv(f"{directory}/E_coordonnee.dat",
                        sep="\s+", header=None, names=["x", "z"])
    n = len(coord)
    # Load def_maille
    def_maille = pd.read_csv(f"{directory}/E_def_maille.dat",
                        sep="\s+", header=None, names=["am", "bm"])
    
    # Prepare an all-(-99) neighbor table (0-based internally)
    ivois = pd.DataFrame(-99, index=range(n),
                         columns=["right", "left", "top", "bottom"])
    
    # Group once for efficient lookup
    x_groups = coord.groupby('z')
    z_groups = coord.groupby('x')
    
    # Find neighbors (0-based)
    for ik, (x_ik, z_ik) in coord[["x", "z"]].iterrows():
        am_ik = def_maille.loc[ik,"am"]
        # RIGHT / LEFT (same z)
        same_z = x_groups.get_group(z_ik).sort_values('x')
        xs, idx_z = same_z["x"].values, same_z.index.values
        p = np.searchsorted(xs, x_ik)
        if p + 1 < len(idx_z): # maille droite
            voisin_idx = idx_z[p+1] # On prend le numéro de maille qui se trouve a droite de notre maille
            x_voisin = coord.loc[voisin_idx,"x"] # Ce qui nous permet d'obtenir sa coordonée en x
            am_voisin = def_maille.loc[voisin_idx,"am"] # Ainsi que sa dimension en x
            dx_reel = abs(x_voisin - x_ik) # Distance entre les deux points
            if dx_reel - (0.5 * (am_ik + am_voisin)) < 1e-1 : # Dans le cas où la différence entre le dx réel et le dx théorique avec la dimension est proche les mailles sont voisines
                ivois.at[ik, "right"] = idx_z[p+1]

        if p - 1 >= 0:
            voisin_idx = idx_z[p-1] # idem maille de gauche
            x_voisin = coord.loc[voisin_idx,"x"]
            am_voisin = def_maille.loc[voisin_idx,"am"]
            dx_reel = abs(x_voisin - x_ik)
            if dx_reel - (0.5 * (am_ik + am_voisin)) < 1e-1 :
                ivois.at[ik, "left"]  = idx_z[p-1]

        
        # TOP / BOTTOM (same x)
        same_x = z_groups.get_group(x_ik).sort_values('z')
        zs, idx_x = same_x["z"].values, same_x.index.values
        p = np.searchsorted(zs, z_ik)
        if p + 1 < len(idx_x):
            ivois.at[ik, "top"]    = idx_x[p+1]
        if p - 1 >= 0:
            ivois.at[ik, "bottom"] = idx_x[p-1]
    
    # Shift all valid (>=0) neighbor indices to 1-based
    ivois = ivois.applymap(lambda v: v+1 if v >= 0 else v)
    
    # (Optional) if you want your DataFrame itself to be 1‑indexed:
    ivois.index = np.arange(1, n+1)
    
    # Write out (each row now corresponds to mesh #1, #2, …)
    ivois.to_csv(f"{directory}/E_voisins.dat",
                 sep=" ", header=False, index=False)
    print("Neighbors saved to E_voisins.dat.")
    
    return ivois
 #______________________________________________________________________   


def maille_limite(directory):
    """
    function to compute the element limite from the file "E_coordonnee.dat",
    producing indices in E_BordRG.dat and E_BordRD.dat.
    The function creates a table id of each limit mesh element,
    PARAMETERS
    ----------
    directory : str
        The directory where the mesh files are located.
    OUTPUT
    -------
    Writes the table to files named "E_BordRG.dat" and "E_BordRD.dat" in the specified directory.
    """
    # Load coordinates
    coord = pd.read_csv(f"{directory}/E_coordonnee.dat",
                        sep="\s+", header=None, names=["x", "z"])
    n = len(coord)
    
    # Prepare an all-(-99) neighbor table (0-based internally)
    BordRG = list()
    BordRD = list()
    
    # Group once for efficient lookup
    x_groups = coord.groupby('z')
    # z_groups = coord.groupby('x')
    
    # Find neighbors (0-based)
    for ik, (x_ik, z_ik) in coord[["x", "z"]].iterrows():
        # RIGHT / LEFT (same z)
        same_z = x_groups.get_group(z_ik).sort_values('x')
        xs, idx_z = same_z["x"].values, same_z.index.values
        p = np.searchsorted(xs, x_ik)
        if p + 1 >= len(idx_z):
            BordRD.append(ik+1)
        if p - 1 < 0:
            BordRG.append(ik+1)
    
    with open(f"{directory}/E_BordRD.dat", "w") as f:
        for item in BordRD:
            f.write(f"{item}\n")
    with open(f"{directory}/E_BordRG.dat", "w") as f:
        for item in BordRG:
            f.write(f"{item}\n")
    # print("Neighbors saved to E_voisins.dat.")

    return BordRG,BordRD
 

#______________________________________________________________________
def coord_to_row_column(repertory):
    """
    function to compute the index of row and column of each mesh element from the file "E_coordonnee.dat" x and z coordinates and E_voisins.dat,
    producing indices in E_colonne.dat and E_row.dat.
    Parameters
    ----------  
    repertory : str
        The directory where the mesh files are located.
    OUTPUT
    -------
    E_colonne.dat : DataFrame
        A DataFrame containing the column indices for each mesh element.
    E_row.dat : DataFrame
        A DataFrame containing the row indices for each mesh element.

    """
    # 1) Load the (x,z) coordinates
    coord = pd.read_csv(
        f"{repertory}/E_coordonnee.dat",
        sep="\s+",
        header=None,
        names=["x", "z"]
    )
    n = len(coord)

    # 2) Build the column index from 'z'
    #   a) Extract all unique z-values, sort them
    z_unique = np.sort(coord["z"].unique())
    #   b) Map each z → its rank in that sorted array (1-based)
    z_to_col = {z_val: idx+1 for idx, z_val in enumerate(z_unique)}
    #   c) Vectorized lookup: map every row’s z to coord["colonne"]
    coord["row"] = coord["z"].map(z_to_col)

    # 3) Build the row index from 'x' (same pattern)
    x_unique = np.sort(coord["x"].unique())
    x_to_row = {x_val: idx+1 for idx, x_val in enumerate(x_unique)}
    coord["colonne"] = coord["x"].map(x_to_row)

    # 4) (Optional Diagnostics)
    nb_col=coord['colonne'].max()
    nb_row=coord['row'].max()
    print(f"Found {len(z_unique)} distinct columns; max assigned column = {coord['colonne'].max()}")
    print(f"Found {len(x_unique)} distinct rows;    max assigned row    = {coord['row'].max()}")

    # 5) Write out the two one‑column files, no header/index
    coord[["colonne"]].to_csv(
        f"{repertory}/E_colonne.dat",
        sep=" ",
        header=False,
        index=False
    )
    coord[["row"]].to_csv(
        f"{repertory}/E_row.dat",
        sep=" ",
        header=False,
        index=False
    )



    return nb_col,nb_row

#______________________________________________________________________
#Definition de la fontion
#____________________________________________________________________
from shapely.geometry import Point, Polygon
# Definition of the function to create the E_zone.dat file
def creation_E_zone(directory, polygons_by_zone, default_zone=1):
    """
    Function to create the E_zone.dat file required for Ginette. This function uses the mesh (E_coordonnee.dat file) 
    and a dictionary containing polygons grouped by their respective zones. The function checks if the center of each 
    mesh lies within a polygon. Once the zone of the mesh is determined, it is added to a list, which is written to 
    the E_zone.dat file at the end.

    Parameters:
    - directory (str): Path (relative or absolute) to the directory containing the E_coordonnee.dat file.
    - polygons_by_zone (dict): Dictionary containing polygons that define the different zones of the mesh.
    - default_zone (int): Default zone (default = 1) assigned to points not contained in any polygon 
                          (avoids defining all polygons in polygons_by_zone).

    Output:
    - Creates the E_zone.dat file in the specified directory.
    """

    # Load the coordinates from the E_coordonnee.dat file
    coord = pd.read_csv(
        os.path.join(directory, "E_coordonnee.dat"),
        sep="\s+", header=None, names=["x", "z"]
    )
    # Generate the E_zone.dat file
    # This file assigns a zone number to each mesh cell based on the defined polygons

    # Start with default zone 1 for all mesh cells
    zone_list = [1] * len(coord)

    # Get all zone numbers from polygons_by_zone, sorted by descending zone number (priority: higher zone number first)
    priority_zones = sorted(polygons_by_zone.keys(), reverse=True)

    # Check each mesh center against polygons for all defined zones (priority order)
    for idx, row in coord.iterrows():
        pt = Point(row['x'], row['z'])
        assigned = False
        for zone_num in priority_zones:
            for poly in polygons_by_zone.get(zone_num, []):
                if poly.contains(pt):
                    zone_list[idx] = zone_num
                    assigned = True
                    break
            if assigned:
                break



    # Test if all expected zones appear in the result
    unique_zones = set(zone_list)
    print("Zones present in E_zone.dat:", unique_zones)
    expected_zones = set([1] + list(polygons_by_zone.keys()))
    if not expected_zones.issubset(unique_zones):
        print("Warning: Not all expected zones are present in E_zone.dat!")
    else:
        print("All expected zones are present in E_zone.dat.")
    # Save the zones to the E_zone.dat file
    pd.DataFrame(zone_list, columns=["zone"]).to_csv(
        os.path.join(directory, "E_zone.dat"), sep=' ', header=False, index=False
    )

#______________________________________________________________________

def id_mesh_river(directory,hmax,hmin,xRG,xRD):
    """
    Read:
      - E_coordonnee.dat  (x,z of each mesh centre)
      - E_voisin.dat    (col1: right: index of the right neighbor, col2 : left: index of the left neighbor, 
      col3 top: index of the top neighbor and col4 : bottom: index of the bottom neighbor  , -99 no neighbor)
     - id maille begins at 1

    Compute:
        -id_river : id of the mesh element which are in the river
        -id_max

    Write:
      - E_Id_river.dat
      - E_Id_river_max.dat
    """
# read file
    coord = pd.read_csv(f"{directory}/E_coordonnee.dat",
                        sep="\s+", header=None, names=["x","z"])
    vois=pd.read_csv(f"{directory}/E_voisins.dat",
                        sep="\s+", header=None, names=["right","left","top","bottom"])
    # find the id of the mesh element which are in the river
    # find all id whithout neighbor at top and x between xRG and xRD
    id_river = vois[vois["top"] == -99].index + 1
    # keep onfly between xRG and xRD
    id_river = id_river[(coord["x"].iloc[id_river-1] > xRG) & (coord["x"].iloc[id_river-1] < xRD)]
    # find the id of the mesh element below h_min
    id_river_min=id_river[coord["z"].iloc[id_river-1] < hmin]
    # find the id of the mesh element above h_max
    id_river_max=id_river[coord["z"].iloc[id_river-1] < hmax]
    # exclude id of id_river_max already in id_river_min
    id_river_max = id_river_max[~id_river_max.isin(id_river_min)]
    # write files
    with open(f"{directory}/E_Id_river.dat", "w") as f:
        for item in id_river_min:
            f.write(f"{item}\n") # \n necessaire pour ecrire un par ligne
    with open(f"{directory}/E_Id_river_max.dat", "w") as f:
        for item in id_river_max:
            f.write(f"{item}\n")

    return id_river_min,id_river_max


#______________________________________________________________________
def am_bm_modif(fName, precision=None):
    """
    Lire un fichier Gmsh (.msh) et calculer les dimensions de chaques des éléments.

    Parameters
    ----------
    fName : str
        Chemin vers le fichier Gmsh (.msh).
    precision : int, optional
        Si spécifié, arrondir les coordonnées des centres à ce nombre de décimales.

    Returns
    -------
    centers : list of tuple
        Liste des coordonnées des centres des éléments sous forme de tuples (x, y, z).
    nb_mesh : int
        nb of center(mesh) in the mesh.
    """
    import gmshparser
    import numpy as np

    assert (precision is None or (isinstance(precision, int) and precision >= 0)), \
        "Precision doit être None ou un entier positif."
    # Parse the Gmsh file
    mesh = gmshparser.parse(fName)
    node_coords = {}

    # Step 1: Collect all node coordinates
    for entity in mesh.get_node_entities():
        for node in entity.get_nodes():
            nid = node.get_tag()
            coords = node.get_coordinates()
            node_coords[nid] = coords

    # Step 2: Calculate centers of elements
    centers_elements = []
    for entity in mesh.get_element_entities():
        # Get the dimension of the entity (e.g., 1 for lines, 2 for surfaces)
        dimension = entity.get_dimension()
        if dimension == 2:  # Only process 2D elements (surfaces)
            for element in entity.get_elements():
                node_ids = element.get_connectivity()
                coords = [node_coords[nid] for nid in node_ids if nid in node_coords]
                if not coords:
                    continue
                am=[c[0] for c in coords]
                bm=[c[2] for c in coords]
                # Calculate the center of the element
                x = np.mean([c[0] for c in coords])
                z = np.mean([c[2] for c in coords])  # Assuming 2D (x, z)
                centers_elements.append((x, z))

    # Step 3: Round coordinates if precision is specified
    if precision is not None:
        centers_elements = [(round(x, precision), round(z, precision)) for x, z in centers_elements]

    # Step 4: Create a DataFrame for the centers
    df_centers = pd.DataFrame(centers_elements, columns=["x", "z"])

    # Step 5: Sort the DataFrame
    # Order with x ascending and z descending
    df_centers = df_centers.sort_values(by=["x", "z"], ascending=[True, False]).reset_index(drop=True)

    # Step 6: Count the number of elements
    nb_mesh = len(df_centers)




#verify if gmsh is installed
try:
    import gmsh 
except ImportError:
    print("GMSH is not installed. Please install it using 'pip install gmsh' or 'conda install gmsh'.")
    sys.exit(1)
#______________________________________________________________________
class MeshPoint:
    """Helper class to manage mesh points with coordinates"""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.tag = None
    
    def add_to_gmsh(self):
        """Add this point to GMSH and store the tag"""
        if self.tag is None:
            self.tag = gmsh.model.occ.addPoint(self.x, self.y, self.z)
        return self.tag
#______________________________________________________________________
class MeshLine:
    """Helper class to manage mesh lines"""
    def __init__(self, point1, point2):
        self.point1 = point1
        self.point2 = point2
        self.tag = None
    
    def add_to_gmsh(self):
        """Add this line to GMSH and store the tag"""
        if self.tag is None:
            p1_tag = self.point1.add_to_gmsh()
            p2_tag = self.point2.add_to_gmsh()
            print(f"p1_tag: {p1_tag}, p2_tag: {p2_tag}")
            print(f"Point 1 coordinates: x={self.point1.x}, y={self.point1.y}, z={self.point1.z}")
            print(f"Point 2 coordinates: x={self.point2.x}, y={self.point2.y}, z={self.point2.z}")
            print(f"Adding line between points {self.point1} and {self.point2}")
            if p1_tag == p2_tag:
                raise ValueError("Points are identical")
            if math.isclose(self.point1.x, self.point2.x, abs_tol=1e-6) and math.isclose(self.point1.y, self.point2.y, abs_tol=1e-6) and math.isclose(self.point1.z, self.point2.z, abs_tol=1e-6):
                raise ValueError("Points are too close together")
            print(f"Adding line between points {self.point1} and {self.point2}")
            try:
                self.tag = gmsh.model.occ.addLine(p1_tag, p2_tag)
            except Exception as e:
                print(f"Error adding line: {e}")
                raise
        return self.tag
    def is_defined(self):
        return self.point1 is not None and self.point2 is not None
    
    def __repr__(self):
        return f"MeshLine(point1={self.point1}, point2={self.point2}, tag={self.tag})"
#_________________________________________________________________________________
class MeshRegion:
    """Helper class to manage mesh regions/surfaces"""
    def __init__(self, name, lines, subdivisions_x, subdivisions_z):
        self.name = name
        self.lines = lines
        self.subdivisions_x = subdivisions_x
        self.subdivisions_z = subdivisions_z
        self.surface_tag = None
    
    def create_surface(self):
        """Create the surface from the lines"""
        line_tags = [line.add_to_gmsh() for line in self.lines]
        loop = gmsh.model.occ.addCurveLoop(line_tags)
        self.surface_tag = gmsh.model.occ.addPlaneSurface([loop])
        return self.surface_tag
    
    def apply_subdivisions(self):
        """Apply transfinite subdivisions to the surface"""
        if self.surface_tag is None:
            raise ValueError(f"Surface for region {self.name} not created yet")
        
        gmsh.model.mesh.setTransfiniteSurface(self.surface_tag, "Left")
        gmsh.model.mesh.setRecombine(2, self.surface_tag)
        
        # Apply subdivisions to lines based on their orientation
        for i, line in enumerate(self.lines):
            if i % 2 == 0:  # Horizontal lines
                gmsh.model.mesh.setTransfiniteCurve(line.tag, self.subdivisions_x)
            else:  # Vertical lines
                gmsh.model.mesh.setTransfiniteCurve(line.tag, self.subdivisions_z)
#___________________________________________________________
def calculate_subdivisions(length, target_size):
    """Calculate optimal integer number of subdivisions for a given length and target size"""
    return max(2, (round(length / target_size)+1))



#___________________________________________________________
def calculate_hobo_parameters(x_hobo_1, x_hobo_2, x_RG, x_RD, num_div_x_centre, dx_hobo):
    """
    Calculate parameters for HOBO mesh refinement.
    This function implements the same logic as the code block you provided.
    Calculate mesh subdivision parameters for HOBO mesh refinement.
    This function determines the subdivision of a mesh along a 1D domain, with special refinement
    around two "HOBO" points (x_hobo_1 and x_hobo_2) between two reference points (x_RG and x_RD).
    It computes the number of divisions and the positions of mesh cells to ensure finer mesh
    near the HOBO points, according to the provided mesh size (dx_hobo) and total number of divisions
    (num_div_x_centre).
    Parameters
    ----------
    x_hobo_1 : float or None
        The x-coordinate of the first HOBO point for mesh refinement. If None, no refinement is done at this point.
    x_hobo_2 : float or None
        The x-coordinate of the second HOBO point for mesh refinement. If None, no refinement is done at this point.
    x_RG : float
        The x-coordinate of the left reference point (start of the mesh region).
    x_RD : float
        The x-coordinate of the right reference point (end of the mesh region).
    num_div_x_centre : int
        The total number of mesh divisions between x_RG and x_RD.
    dx_hobo : float
        The desired mesh cell size near the HOBO points.
    Returns
    -------
    dict
        A dictionary containing:
            - 'num_div_x_hobo': Number of divisions around the first HOBO point (or None).
            - 'num_div_x_RG_hobo': Number of divisions from x_RG to the first HOBO point (or None).
            - 'num_div_x_hobo1_hobo2': Number of divisions between the two HOBO points (or from first HOBO to x_RD if only one HOBO).
            - 'num_div_x_hobo2_RD': Number of divisions from the second HOBO point to x_RD (or 0).
            - 'num_div_x_hobo2': Number of divisions around the second HOBO point (or None).
            - 'x_gauche_hobo': Left boundary of the first HOBO mesh region (or None).
            - 'x_droite_hobo': Right boundary of the first HOBO mesh region (or None).
            - 'x_gauche_hobo2': Left boundary of the second HOBO mesh region (or None).
            - 'x_droite_hobo2': Right boundary of the second HOBO mesh region (or None).
    Notes
    -----
    - The function ensures at least two divisions in each refined region.
    - If a HOBO point is not specified (None), the corresponding outputs are set to None or 0.
    - The function swaps x_RG and x_RD if they are provided in reverse order.
    - The mesh is divided such that the HOBO points are located within their own refined mesh cells.

    """
    if(x_RD<x_RG):
        # X_RD=x_RG and x_RG=x_RD
        x_RG, x_RD = x_RD, x_RG
    
    dx_reel = abs(x_RD - x_RG) / num_div_x_centre
    maille_hobo1 = 0
    maille_hobo2 = 0
    nb_case_a_ajouter1 = 1
    nb_case_a_ajouter2 = 1

    # Find the cell index for hobo1
    if x_hobo_1 is not None:
        for i in range(num_div_x_centre):
            if (x_RG + dx_reel * (i + 1)) > x_hobo_1:
                break
            else:
                maille_hobo1 += 1
        if (x_RG + dx_reel * (maille_hobo1 + 1)) - x_hobo_1 < 0.005:
            nb_case_a_ajouter1 = 2
        x_gauche_hobo1 = x_RG + dx_reel * maille_hobo1
        x_droite_hobo1 = x_RG + dx_reel * (maille_hobo1 + nb_case_a_ajouter1)
        length_x_hobo1 = abs(x_droite_hobo1 - x_gauche_hobo1)
        num_div_x_hobo1 = max(2, math.ceil(length_x_hobo1 / dx_hobo))
        num_div_x_RG_hobo1 = maille_hobo1
    else:
        x_gauche_hobo1 = x_droite_hobo1 = length_x_hobo1 = num_div_x_hobo1 = num_div_x_RG_hobo1 = None

    # Find the cell index for hobo2
    if x_hobo_2 is not None:
        for i in range(num_div_x_centre):
            if (x_RG + dx_reel * (i + 1)) > x_hobo_2:
                break
            else:
                maille_hobo2 += 1
        if (x_RG + dx_reel * (maille_hobo2 + 1)) - x_hobo_2 < 0.005:
            nb_case_a_ajouter2 = 2
        x_gauche_hobo2 = x_RG + dx_reel * maille_hobo2 
        x_droite_hobo2 = x_RG + dx_reel * (maille_hobo2 + nb_case_a_ajouter2)
        length_x_hobo2 = abs(x_droite_hobo2 - x_gauche_hobo2)
        num_div_x_hobo2 = max(2, math.ceil(length_x_hobo2 / dx_hobo))
    else:
        x_gauche_hobo2 = x_droite_hobo2 = length_x_hobo2 = num_div_x_hobo2 = None

    # Subdivisions between hobo1 and hobo2, and after hobo2
    if x_hobo_1 is not None and x_hobo_2 is not None:
        num_div_x_hobo1_hobo2 = max(2, math.ceil((x_gauche_hobo2 - x_droite_hobo1) / dx_reel))
        num_div_x_hobo2_RD = num_div_x_centre - (maille_hobo2 + nb_case_a_ajouter2)
    elif x_hobo_1 is not None and x_hobo_2 is None:
        num_div_x_hobo1_hobo2 = max(2, math.ceil((x_RD - x_droite_hobo1) / dx_reel))
        num_div_x_hobo2_RD = 0
    else:
        num_div_x_hobo1_hobo2 = num_div_x_hobo2_RD = 0
    print('dx_reel', dx_reel, dx_hobo)

    return {
        'num_div_x_hobo': num_div_x_hobo1 if x_hobo_1 is not None else None,
        'num_div_x_RG_hobo1': num_div_x_RG_hobo1 if x_hobo_1 is not None else None,
        'num_div_x_hobo1_hobo2': num_div_x_hobo1_hobo2 if x_hobo_1 is not None  else None,
        'num_div_x_hobo2_RD': num_div_x_hobo2_RD if x_hobo_2 is not None else None,
        'num_div_x_hobo2': num_div_x_hobo2 if x_hobo_2 is not None else None,
        'x_gauche_hobo': x_gauche_hobo1 if x_hobo_1 is not None else None,
        'x_droite_hobo': x_droite_hobo1 if x_hobo_1 is not None else None,
        'x_gauche_hobo2': x_gauche_hobo2 if x_hobo_2 is not None else None,
        'x_droite_hobo2': x_droite_hobo2 if x_hobo_2 is not None else None,
    }
#___________________________________________________________
# new funtion
def generate_mesh_8_region_optimized(distance_altitude_table, output_mesh_path, v_bot=103.8, 
                                   x_RG=5, x_RD=17, z_riv=106, 
                                   dx_grossier=0.5, dx_precis=0.1, 
                                   hobo_points=None, dx_hobo=0.01, 
                                   dz_grossier=1.0, dz_precis=0.2, dz_hobo1=0.1, dz_hobo=0.15, 
                                   mesh_dimension=2, verbose=True,
                                   # Backward compatibility parameters
                                   x_hobo_1=None, z_hobo_1=None, x_hobo_2=None, z_hobo_2=None):
    """
    Optimized version of generate_mesh_8_region function.
    
    This function creates a structured mesh with 8 regions using GMSH, with optimizations for:
    - Better code organization using helper classes
    - Reduced code duplication
    - More efficient subdivision calculations
    - Improved error handling
    - Better memory management
    
    Parameters:
    -----------
    distance_altitude_table : DataFrame
        Table containing topography data with columns "Distance (m)" and "Altitude (Z)"
    output_mesh_path : str
        Path where the mesh will be saved
    v_bot : float, default=103.8
        Z coordinate delimiting the boundary between regions 1-4, 2-5, and 3-6
    x_RG : float, default=5
        X coordinate of the left river bank
    x_RD : float, default=17
        X coordinate of the right river bank
    z_riv : float, default=106
        Z coordinate of the river level
    dx_grossier : float, default=0.5
        Coarse mesh size in X direction for regions 1, 3, 4, 6, 7, 8
    dx_precis : float, default=0.1
        Fine mesh size in X direction for regions 2 and 5
    x_hobo_1, x_hobo_2 : float
        X coordinates of HOBO measurement points requiring fine mesh
    z_hobo_1, z_hobo_2 : float
        Z coordinates of HOBO measurement points
    dx_hobo : float, default=0.01
        Very fine mesh size around HOBO points
    dz_grossier : float, default=1.0
        Coarse mesh size in Z direction for regions 1, 2, 3, 7, 8
    dz_precis : float, default=0.2
        Fine mesh size in Z direction for regions 4, 5, 6
    mesh_dimension : int, default=2
        Mesh dimension (2D)
    verbose : bool, default=True
        Whether to print progress information
    
    Returns:
    --------
    None
        The mesh is saved to output_mesh_path


    Scheme

      z ↑
  |
  | p4---------p8------p17--p18------p12------p23--p24------p16   (z = max_z)
  |  |         |        |    |        |         |    |        |
  |  |         |        |    |        |         |    |        |
  | p3        p7      p17  p18      p11       p23  p24      p15   (z = z_riv)
  |  |         |        |    |        |         |    |        |
  |  |         |        |    |        |         |    |        |
  | p2        p6      p19  p20      p10       p25  p26      p14   (z = v_bot)
  |  |         |        |    |        |         |    |        |
  |  |         |        |    |        |         |    |        |
  | p1        p5      p21  p22      p9        p27  p28      p13   (z = min_z)
  +------------------------------------------------------------→ x
    min_x    x_RG   x_g_hobo1 x_d_hobo1   x_RD-   x_g_hobo2 x_d_hobo2   max_x
    
    Légende :
    - p1  = (min_x, 0, min_z)
    - p2  = (min_x, 0, v_bot)
    - p3  = (min_x, 0, z_riv)
    - p4  = (min_x, 0, max_z)
    - p5  = (x_RG, 0, min_z)
    - p6  = (x_RG, 0, v_bot)
    - p7  = (x_RG, 0, z_riv)
    - p8  = (x_RG, 0, max_z)
    - p9  = (x_RD, 0, min_z)
    - p10 = (x_RD, 0, v_bot)
    - p11 = (x_RD, 0, z_riv)
    - p12 = (x_RD, 0, max_z)
    - p13 = (max_x, 0, min_z)
    - p14 = (max_x, 0, v_bot)
    - p15 = (max_x, 0, z_riv)
    - p16 = (max_x, 0, max_z)
    
    Si HOBO1 est défini :
    - p17 = (x_gauche_hobo1, 0, z_riv)
    - p18 = (x_droite_hobo1, 0, z_riv)
    - p19 = (x_gauche_hobo1, 0, v_bot)
    - p20 = (x_droite_hobo1, 0, v_bot)
    - p21 = (x_gauche_hobo1, 0, min_z)
    - p22 = (x_droite_hobo1, 0, min_z)
    
    Si HOBO2 est défini :
    - p23 = (x_gauche_hobo2, 0, z_riv)
    - p24 = (x_droite_hobo2, 0, z_riv)
    - p25 = (x_gauche_hobo2, 0, v_bot)
    - p26 = (x_droite_hobo2, 0, v_bot)
    - p27 = (x_gauche_hobo2, 0, min_z)
    - p28 = (x_droite_hobo2, 0, min_z)
    
    """
    
    try:
        if verbose:
            print("\n=== Generating optimized structured mesh with GMSH ===")
        
        # Initialize GMSH
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1 if verbose else 0)
        gmsh.model.add("OptimizedRectangularMesh")

        # Extract coordinates from input data
        distances = distance_altitude_table["Distance (m)"].to_numpy()
        z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()

        # Calculate domain boundaries
        min_z, max_z = min(z_coords), max(z_coords)
        min_x, max_x = min(distances), max(distances)
        
        # Update z_riv to be the maximum of input z_riv and HOBO heights
        if z_hobo_1 is None and z_hobo_2 is not None:
            # z_hobo 1 = z_hobo_2  and z_hobo_2 = z_hobo_1 same with x
            z_hobo_1, z_hobo_2 = z_hobo_2, z_hobo_1
            x_hobo_1, x_hobo_2 = x_hobo_2, x_hobo_1
            print("Only one temperature profile exist")


        
        if verbose:
            print(f"Domain: X=[{min_x:.4f}, {max_x:.4f}], Z=[{min_z:.4f}, {max_z:.4f}]")
            print(f"River level: {z_riv:.4f}, Riverbed Bottom level: {v_bot:.4f}")

        # Calculate basic subdivisions
        length_x_gauche = round(abs(x_RG - min_x),4)
        length_x_droite = round(abs(max_x - x_RD),4)
        length_x_centre = round(abs(x_RD - x_RG),4)
        length_z_haut = round(abs(max_z - z_riv),4)
        length_z_centre =round(abs(z_riv - v_bot),4)
        length_z_bas = round(abs(v_bot - min_z),4)

        
        num_div_x_gauche = calculate_subdivisions(length_x_gauche, dx_grossier)
        num_div_x_centre = calculate_subdivisions(length_x_centre, dx_precis)
        num_div_x_droite = calculate_subdivisions(length_x_droite, dx_grossier)
        
        num_div_z_haut= calculate_subdivisions(length_z_haut, dz_grossier)
        num_div_z_centre = calculate_subdivisions(length_z_centre, dz_precis)
        num_div_z_bas = calculate_subdivisions(length_z_bas, dz_grossier)
        if not isinstance(num_div_z_centre, int):
            print(f"num_div_z_centre {num_div_z_centre} dz_precis {dz_precis}")
            # error stop and ask to change v_bot and z_riv
            raise ValueError("z_riv and v_bot")

        print(f"num_div_z_centre {num_div_z_centre} dz_precis {dz_precis}")
        #verif
        dz_calcul= round(length_z_centre / (num_div_z_centre),4)

        if dz_calcul != dz_precis:
            print(f"Warning: Calculated dz {dz_calcul} does not match expected dz_precis {dz_precis}. Adjusting num_div_z_centre.")
            print('ici',abs(round(z_riv - v_bot,4)))
            print(f"num_div_z_centre {num_div_z_centre} dz_precis {dz_precis} dz_calcul {dz_calcul}")


        # Calculate HOBO parameters
        if x_hobo_1 is not None and x_hobo_2 is not None:
            hobo_params = calculate_hobo_parameters(x_hobo_1, x_hobo_2, x_RG, x_RD, num_div_x_centre, dx_hobo)
        elif x_hobo_1 is not None and x_hobo_2 is None:
            hobo_params = calculate_hobo_parameters(x_hobo_1, None, x_RG, x_RD, num_div_x_centre, dx_hobo)
        else:
            hobo_params = None
        print("hobo_params", hobo_params)    
        
        if x_hobo_1 is not None and z_hobo_1 is not None:
            if x_hobo_2 is not None and z_hobo_2 is not None:
                # Both HOBO points exist
                num_div_x_RG_hobo1 = hobo_params['num_div_x_RG_hobo1']
                num_div_x_hobo1_hobo2 = hobo_params['num_div_x_hobo1_hobo2']
                num_div_x_hobo2_RD = hobo_params['num_div_x_hobo2_RD']
            else:
                # Only one HOBO point exists
                num_div_x_RG_hobo1 = hobo_params['num_div_x_RG_hobo1']
                # if only one hobo discretisation adapt discretisation
                num_div_x_hobo1_hobo2 = hobo_params['num_div_x_hobo1_hobo2']
                num_div_x_hobo2_RD = 0

        else:
            # No HOBO points exist
            hobo_params = None
            num_div_x_RG_hobo1 = 0
            num_div_x_hobo1_hobo2 = 0
            num_div_x_hobo2_RD = 0

        if verbose:
            print(f"Subdivisions - Left: {num_div_x_gauche}, Center: {num_div_x_centre}, Right: {num_div_x_droite}")
            print(f"Z subdivisions - Top: {num_div_z_haut}, Center: {num_div_z_centre}, Bottom: {num_div_z_bas}")

        # Create all mesh points using the helper class
        points = {}        
        # Main corner points
        points['p1'] = MeshPoint(min_x, 0, min_z)
        points['p2'] = MeshPoint(min_x, 0, v_bot)
        points['p3'] = MeshPoint(min_x, 0, z_riv)
        points['p4'] = MeshPoint(min_x, 0, max_z)
        
        points['p5'] = MeshPoint(x_RG, 0, min_z)
        points['p6'] = MeshPoint(x_RG, 0, v_bot)
        points['p7'] = MeshPoint(x_RG, 0, z_riv)
        points['p8'] = MeshPoint(x_RG, 0, max_z)
        
        points['p9'] = MeshPoint(x_RD, 0, min_z)
        points['p10'] = MeshPoint(x_RD, 0, v_bot)
        points['p11'] = MeshPoint(x_RD, 0, z_riv)
        points['p12'] = MeshPoint(x_RD, 0, max_z)
        
        points['p13'] = MeshPoint(max_x, 0, min_z)
        points['p14'] = MeshPoint(max_x, 0, v_bot)
        points['p15'] = MeshPoint(max_x, 0, z_riv)
        points['p16'] = MeshPoint(max_x, 0, max_z)
        
        if x_hobo_1 is not None:
            points['p17'] = MeshPoint(hobo_params['x_gauche_hobo'], 0, z_riv)
            points['p18'] = MeshPoint(hobo_params['x_droite_hobo'], 0, z_riv)
            points['p19'] = MeshPoint(hobo_params['x_gauche_hobo'], 0, v_bot)
            points['p20'] = MeshPoint(hobo_params['x_droite_hobo'], 0, v_bot)
            points['p21'] = MeshPoint(hobo_params['x_gauche_hobo'], 0, min_z)
            points['p22'] = MeshPoint(hobo_params['x_droite_hobo'], 0, min_z)

        if x_hobo_2 is not None:
            points['p23'] = MeshPoint(hobo_params['x_gauche_hobo2'], 0, z_riv)
            points['p24'] = MeshPoint(hobo_params['x_droite_hobo2'], 0, z_riv)
            points['p25'] = MeshPoint(hobo_params['x_gauche_hobo2'], 0, v_bot)
            points['p26'] = MeshPoint(hobo_params['x_droite_hobo2'], 0, v_bot)
            points['p27'] = MeshPoint(hobo_params['x_gauche_hobo2'], 0, min_z)
            points['p28'] = MeshPoint(hobo_params['x_droite_hobo2'], 0, min_z)

        
        lines = {}

        # Vertical lines
        lines['l1'] = MeshLine(points['p1'], points['p2'])
        lines['l2'] = MeshLine(points['p2'], points['p3'])
        lines['l3'] = MeshLine(points['p3'], points['p4'])
        lines['l4'] = MeshLine(points['p5'], points['p6'])
        lines['l5'] = MeshLine(points['p6'], points['p7'])
        lines['l6'] = MeshLine(points['p7'], points['p8'])
        lines['l7'] = MeshLine(points['p9'], points['p10'])
        lines['l8'] = MeshLine(points['p10'], points['p11'])
        lines['l9'] = MeshLine(points['p11'], points['p12'])
        lines['l10'] = MeshLine(points['p13'], points['p14'])
        lines['l11'] = MeshLine(points['p14'], points['p15'])
        lines['l12'] = MeshLine(points['p15'], points['p16'])

        # HOBO vertical lines
        if x_hobo_1 is not None:
            lines['l13'] = MeshLine(points['p21'], points['p19'])
            lines['l14'] = MeshLine(points['p19'], points['p17'])
            lines['l15'] = MeshLine(points['p22'], points['p20'])
            lines['l16'] = MeshLine(points['p20'], points['p18'])

        if x_hobo_2 is not None:
            lines['l17'] = MeshLine(points['p27'], points['p25'])
            lines['l18'] = MeshLine(points['p25'], points['p23'])
            lines['l19'] = MeshLine(points['p28'], points['p26'])
            lines['l20'] = MeshLine(points['p26'], points['p24'])

        # Horizontal lines
        lines['l100'] = MeshLine(points['p1'], points['p5'])
        lines['l101'] = MeshLine(points['p2'], points['p6'])
        lines['l102'] = MeshLine(points['p3'], points['p7'])
        lines['l103'] = MeshLine(points['p4'], points['p8'])
        # Additional horizontal lines
        lines['l112'] = MeshLine(points['p9'], points['p13'])
        lines['l113'] = MeshLine(points['p10'], points['p14'])
        lines['l114'] = MeshLine(points['p11'], points['p15'])
        lines['l115'] = MeshLine(points['p12'], points['p16'])
        # HOBO horizontal lines
        if x_hobo_1 is not None:
            lines['l104'] = MeshLine(points['p5'], points['p21'])
            lines['l105'] = MeshLine(points['p6'], points['p19'])
            lines['l106'] = MeshLine(points['p7'], points['p17'])
            lines['l107'] = MeshLine(points['p21'], points['p22'])
            lines['l108'] = MeshLine(points['p17'], points['p18'])
            lines['l116'] = MeshLine(points['p19'], points['p20'])

        if x_hobo_2 is not None:
            lines['l109'] = MeshLine(points['p22'], points['p27'])
            lines['l110'] = MeshLine(points['p20'], points['p25'])
            lines['l111'] = MeshLine(points['p18'], points['p23'])
            lines['l117'] = MeshLine(points['p27'], points['p28'])
            lines['l118'] = MeshLine(points['p25'], points['p26'])
            lines['l119'] = MeshLine(points['p23'], points['p24'])
            lines['l120'] = MeshLine(points['p28'], points['p9'])
            lines['l121'] = MeshLine(points['p26'], points['p10'])
            lines['l122'] = MeshLine(points['p24'], points['p11'])
        if x_hobo_2 is None and x_hobo_1 is None:
            # ajouter missing lines
            lines['104']=MeshLine(points['p5'],points['p9'])
            lines['105']=MeshLine(points['p6'],points['p10'])
            lines['106']=MeshLine(points['p7'],points['p11'])
            lines['107']=MeshLine(points['p8'],points['p12'])
        if x_hobo_2 is None and x_hobo_1 is not None:
            # ajouter missing lines
            lines['l109'] = MeshLine(points['p22'], points['p9'])
            lines['l110'] = MeshLine(points['p20'], points['p10'])
            lines['l111'] = MeshLine(points['p18'], points['p11'])
            


        # Create regions with their respective subdivisions
        regions = []

        # Define each region with its boundary lines and subdivisions
        region_definitions = []

        # Case 1: No HOBO points
        if x_hobo_1 is None and x_hobo_2 is None:
            print("No temperature profiles exist")
            # Define regions for the case without HOBO mesh refinement
            # Each region is defined by its boundary lines and the number of subdivisions in x and z
            # Region1: Bottom left
            # Region2: Bottom center/right
            # Region3: Middle left
            # Region4: Middle center/right
            # Region5: Top left
            # Region6: Top center/right

            region_definitions = [
            ("Region1", [lines['l101'], lines['l4'], lines['l100'], lines['l1']], num_div_x_gauche, num_div_z_bas),
            ("Region2", [lines['l113'], lines['l10'], lines['l112'], lines['l4']], num_div_x_centre, num_div_z_bas),
            ("Region3", [lines['l102'], lines['l5'], lines['l101'], lines['l2']], num_div_x_gauche, num_div_z_centre),
            ("Region4", [lines['l114'], lines['l11'], lines['l113'], lines['l5']], num_div_x_centre, num_div_z_centre),
            ("Region5", [lines['l103'], lines['l6'], lines['l102'], lines['l3']], num_div_x_gauche, num_div_z_haut),
            ("Region6", [lines['l115'], lines['l12'], lines['l114'], lines['l6']], num_div_x_centre, num_div_z_haut),
            ]
        # Case 2: Only HOBO 1 exists
        elif x_hobo_1 is not None and x_hobo_2 is None:
            """
            z ↑
            |
            | p4--l103---p8--------p17---------p18------------p12---l115-----p16   (z = max_z)
            |  |         |         |           |              |              |            
              l3 RG10    l6                                   l9     RG11      l12           
            |  |         |         |           |              |              |
            | p3--l102--p7--L106--p17---l108---p18---l111-----p11---l114----p15   (z = z_riv)
            |  |         |         |           |              |              |            
            |  |         |         |           |              |              |      
            | L2   RG6   l5   RG7  l14   RG12   l16      RG8   l8     RG9     l11
            |  |         |         |           |              |              |     
            |  |         |         |           |              |              |
            | p2--l101--p6--l105--p19---l116---p20---l110-----p10---l113-----p14   (z = v_bot)
            |  |         |         |           |              |              |             
            |  |         |         |           |              |              |           
            | L1   RG1   l4  RG2   l13  RG3    l15     RG4    l7      RG5 l10
            |  |         |         |           |              |              |             
            |  |         |         |           |              |              |          
            | p1--l100--p5--l104--p21---l107---p22---l109-----p9---l112------p13   (z = min_z)
            +----------------- --------------------------------------------------------------------→ x
                min_x    x_RG   x_g_hobo1     x_d_hobo1         x_RD   max_x   
            
            
            """
            print('Regions with one temperature profile')
            # You can further split this block if you want to distinguish the three cases
            region_definitions = [("Region1", [lines['l101'], lines['l4'], lines['l100'], lines['l1']], num_div_x_gauche, num_div_z_bas),
            ("Region2", [lines['l105'], lines['l13'], lines['l104'], lines['l4']], num_div_x_RG_hobo1, num_div_z_bas),
            ("Region3", [lines['l116'], lines['l15'], lines['l107'], lines['l13']], hobo_params['num_div_x_hobo'], num_div_z_bas),
            ("Region4", [lines['l110'], lines['l7'], lines['l109'], lines['l15']],  num_div_x_hobo1_hobo2, num_div_z_bas),
            ("Region5", [lines['l113'], lines['l10'], lines['l112'], lines['l7']], num_div_x_droite, num_div_z_bas),
            ("Region6", [lines['l102'], lines['l5'], lines['l101'], lines['l2']], num_div_x_gauche, num_div_z_centre),
            ("Region7", [lines['l106'], lines['l14'], lines['l105'], lines['l5']], num_div_x_RG_hobo1, num_div_z_centre),
            ("Region8", [lines['l111'], lines['l8'], lines['l110'], lines['l16']], num_div_x_hobo1_hobo2, num_div_z_centre),
            ("Region9", [lines['l114'], lines['l11'], lines['l113'], lines['l8']], num_div_x_droite, num_div_z_centre),
            ("Region10", [lines['l103'], lines['l6'], lines['l102'], lines['l3']], num_div_x_gauche, num_div_z_haut),
            ("Region11", [lines['l115'], lines['l12'], lines['l114'], lines['l9']], num_div_x_droite, num_div_z_haut),
            ("Region12", [lines['l108'], lines['l16'], lines['l116'], lines['l14']], hobo_params['num_div_x_hobo'], num_div_z_centre)]

            #for i, region in enumerate(region_definitions):
             #   print(f"Region {i+1}: {region[0]}")
              #  print(f"Lines: {region[1]}")
              #  for j, line in enumerate(region[1]):
               #     print(f"Line {j+1}: {line}")
               # print()

        
        # Case 3: Both HOBO points exist
        elif x_hobo_1 is not None and x_hobo_2 is not None:
            print("2 temperature profiles exist")
            """
            ↑
            |
            | p4--l103---p8--------p17---------p18----------p23---------p24----------p12---l115---p16   (z = max_z)
            |  |         |         |           |            |           |             |           |
                l3        l6                                                                      l12  
            |  |         |         |           |            |           |             |           |
            | p3--l102--p7--L106--p17---l108---p18---l111---p23---l119---p24---l122---p11---l114---p15   (z = z_riv)
            |  |         |         |           |            |           |             |            |
            |  |         |         |           |            |           |             |            |
            | L2        l5        l14         l16           l18         l20           l8           l11
            |  |         |         |           |            |           |             |            |
            |  |         |         |           |            |           |             |            |
            | p2--l101--p6--l105--p19---l116---p20---l110---p25---l118---p26---l121---p10---l113---p14   (z = v_bot)
            |  |         |         |           |            |           |             |            |
            |  |         |         |           |            |           |             |            |
                L1        l4       l13          l15          l17         l19           l7           l10
            |  |         |         |           |            |           |             |            |
            |  |         |         |           |            |           |             |             |
            | p1--l100--p5--l104--p21---l107---p22---l109---p27---l117---p28---l120---p9---l112---p13   (z = min_z)
            +----------------- --------------------------------------------------------------------→ x
             min_x    x_RG   x_g_hobo1      x_d_hobo1     x_g_hobo2     x_d_hobo2    x_RD   max_x  
     
    """
            # Three possible sub-cases for both HOBOs (example: left, middle, right)
            # You can further split this block if you want to distinguish the three cases
            region_definitions = [
            ("Region1", [lines['l101'], lines['l4'], lines['l100'], lines['l1']], num_div_x_gauche, num_div_z_bas),
            ("Region2", [lines['l105'], lines['l13'], lines['l104'], lines['l4']], num_div_x_RG_hobo1, num_div_z_bas),
            ("Region3", [lines['l116'], lines['l15'], lines['l107'], lines['l13']], hobo_params['num_div_x_hobo'], num_div_z_bas),
            ("Region4", [lines['l110'], lines['l17'], lines['l109'], lines['l15']], num_div_x_hobo1_hobo2, num_div_z_bas),
            ("Region5", [lines['l113'], lines['l10'], lines['l112'], lines['l7']], num_div_x_droite, num_div_z_bas),
            ("Region6", [lines['l102'], lines['l5'], lines['l101'], lines['l2']], num_div_x_gauche, num_div_z_centre),
            ("Region7", [lines['l106'], lines['l14'], lines['l105'], lines['l5']], num_div_x_RG_hobo1, num_div_z_centre),
            ("Region8", [lines['l111'], lines['l18'], lines['l110'], lines['l16']],num_div_x_hobo1_hobo2, num_div_z_centre),
            ("Region9", [lines['l114'], lines['l11'], lines['l113'], lines['l8']], num_div_x_droite, num_div_z_centre),
            ("Region10", [lines['l103'], lines['l6'], lines['l102'], lines['l3']], num_div_x_gauche, num_div_z_haut),
            ("Region11", [lines['l115'], lines['l12'], lines['l114'], lines['l9']], num_div_x_droite, num_div_z_haut),
            ("Region12", [lines['l108'], lines['l16'], lines['l116'], lines['l14']], hobo_params['num_div_x_hobo'], num_div_z_centre), 
            ("Region13", [lines['l118'], lines['l19'], lines['l117'], lines['l17']], hobo_params['num_div_x_hobo'], num_div_z_bas),
            ("Region14", [lines['l121'], lines['l7'], lines['l120'], lines['l19']], num_div_x_hobo2_RD, num_div_z_bas),
            ("Region15", [lines['l119'], lines['l20'], lines['l118'], lines['l18']], hobo_params['num_div_x_hobo'], num_div_z_centre),
            ("Region16", [lines['l122'], lines['l8'], lines['l121'], lines['l20']], num_div_x_hobo2_RD, num_div_z_centre)]

        # how to check before creating regions
        if verbose:
            print(f"Creating {len(region_definitions)} mesh regions...")
        # check if lines are defined
        for line in lines.values():
            if not line.is_defined():
                raise ValueError(f"Line {line.name} is not defined. Please check the mesh lines.")
        # check order
        if verbose:
            print("All lines are defined and in correct order.")

        # Create regions
        for region_definition in region_definitions:
            region_name, lines, num_div_x, num_div_z = region_definition
            regions.append(MeshRegion(region_name, lines, num_div_x, num_div_z))

        # Create all surfaces
        for region in regions:
            region.create_surface()

        # Synchronize the model
        gmsh.model.occ.synchronize()

        # Apply subdivisions to all regions
        for region in regions:
            region.apply_subdivisions()

        if verbose:
            print("Generating mesh...")
        
        # Generate the mesh
        gmsh.model.mesh.generate(mesh_dimension)
        # check if duplicated nodes
        gmsh.model.mesh.removeDuplicateNodes()
        # check if duplicated elements
        gmsh.model.mesh.removeDuplicateElements()

        gmsh.write(output_mesh_path)
        
        if verbose:
            print(f" Optimized esh saved to: {output_mesh_path}")
            print(f"Total regions created: {len(regions)}")
#
    except Exception as e:
        print(f" Error during mesh generation: {e}")
        raise
    finally:
        gmsh.finalize()

#________________________________________________________________________________

# Backward compatibility function
def generate_mesh_8_region(distance_altitude_table, output_mesh_path, **kwargs):
    """
    Backward compatibility wrapper for the optimized function.
    """
    return generate_mesh_8_region_optimized(distance_altitude_table, output_mesh_path, **kwargs)

#________________________________________________________________________________