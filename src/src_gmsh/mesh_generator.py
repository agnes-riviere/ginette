import os
import gmsh
import numpy as np
import pandas as pd
import pyvista as pv
import math
import polytools
import quality
import grid
import mapping


# Configuration du répertoire de travail
script_dir = os.path.dirname(os.path.abspath(__file__))  # Répertoire du script
sim_dir = os.path.join(script_dir, "meshdata")  # Chemin absolu vers "meshdata"
Station = "IntD"
altitude_min = 103

if not os.path.exists(sim_dir):
    os.makedirs(sim_dir)

# Chemin vers le fichier de nivellement
file_path = f"/home/ariviere/Documents/Bassin_Orgeval/Donnee_Orgeval_Mines/raw_data/DESC_data/DATA_STATION/{Station}/{Station}.csv"

# Vérification de l'existence du fichier
if not os.path.exists(file_path):
    raise FileNotFoundError(f"Le fichier {file_path} est introuvable !")

# Lecture des données
try:
    data = pd.read_csv(file_path, delim_whitespace=True, header=None)
    x_coords, y_coords, z_coords = data[0].to_numpy(), data[1].to_numpy(), data[2].to_numpy()
except Exception as e:
    raise ValueError(f"Erreur lors de la lecture du fichier : {e}")

# Calculer les distances cumulées
distances = np.zeros(len(x_coords))
for i in range(1, len(x_coords)):
    distances[i] = distances[i - 1] + np.sqrt((x_coords[i] - x_coords[i - 1])**2 + (y_coords[i] - y_coords[i - 1])**2)

distances = np.round(distances, 1)

# Ajout des points aux extrémités
distances = np.insert(distances, 0, 0)
z_coords = np.insert(z_coords, 0, altitude_min)
distances = np.append(distances, distances[-1])
z_coords = np.append(z_coords, altitude_min)
z_coords = np.round(z_coords, 1)

# Création du DataFrame
distance_altitude_table = pd.DataFrame({"Distance (m)": distances, "Altitude (Z)": z_coords})
# plot des données
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.plot(distance_altitude_table["Distance (m)"], distance_altitude_table["Altitude (Z)"], marker='o', linestyle='-')
plt.title("Distance vs Altitude")
plt.xlabel("Distance (m)")
plt.ylabel("Altitude (Z)")
plt.grid()
plt.show()

# Sauvegarde en CSV
output_csv_path = os.path.join(sim_dir, "distance_vs_altitude.csv")
distance_altitude_table.to_csv(output_csv_path, index=False)
def generate_quadrilateral_mesh(distance_altitude_table, output_mesh_path, dx=0.1, dz=0.1, mesh_dimension=2):
    try:
        print("\nGénération du maillage structuré avec Gmsh...")
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("RectangularMesh")

        # Extraction des distances et altitudes
        distances = distance_altitude_table["Distance (m)"].to_numpy()
        z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()
        print(f"Distances : {distances}")
        print(f"Altitudes : {z_coords}")
        # Points du rectangle
        p1 = gmsh.model.occ.addPoint(distances[0], 0, z_coords[0])
        p2 = gmsh.model.occ.addPoint(distances[-1], 0, z_coords[-1])
        p3 = gmsh.model.occ.addPoint(distances[-2], 0, z_coords[-2])
        p4 = gmsh.model.occ.addPoint(distances[1], 0, z_coords[1])

        # Création des lignes et de la surface
        l1 = gmsh.model.occ.addLine(p1, p2)
        l2 = gmsh.model.occ.addLine(p2, p3)
        l3 = gmsh.model.occ.addLine(p3, p4)
        l4 = gmsh.model.occ.addLine(p4, p1)
        loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
        surface = gmsh.model.occ.addPlaneSurface([loop])

        # Synchronisation du modèle
        gmsh.model.occ.synchronize()

        # Calcul des subdivisions
        length_x = distances[-1] - distances[0]
        length_z = max(z_coords) - min(z_coords)
        num_div_x = max(2, math.ceil(length_x / dx))
        num_div_z = max(2, math.ceil(length_z / dz))
        print(f"Nombre de divisions en X : {num_div_x}, Nombre de divisions en Z : {num_div_z}")

        # Application des subdivisions et recombinaison pour quadrilatères
        gmsh.model.mesh.setTransfiniteCurve(l1, num_div_x)
        gmsh.model.mesh.setTransfiniteCurve(l3, num_div_x)
        gmsh.model.mesh.setTransfiniteCurve(l2, num_div_z)
        gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z)
        gmsh.model.mesh.setTransfiniteSurface(surface, "Left")
        gmsh.model.mesh.setRecombine(2, surface)

        # Génération du maillage
        gmsh.model.mesh.generate(mesh_dimension)
        gmsh.write(output_mesh_path)
        print(f"✅ Maillage sauvegardé sous : {output_mesh_path}")

    except Exception as e:
        print(f"Erreur lors de la génération du maillage : {e}")
    finally:
        gmsh.finalize()

# Appel de la fonction
output_mesh_path = os.path.join(sim_dir, "rectangular_mesh.msh")
#generate_quadrilateral_mesh(distance_altitude_table, output_mesh_path, dx=0.1, dz=0.5, mesh_dimension=2)



def generate_quadrilateral_mesh_with_aligned_lines(output_mesh_path, dx_outer=1.0, dz_outer=1.0, dx_inner=0.2, dz_inner=0.2, mesh_dimension=2):
    try:
        print("\nGénération du maillage structuré avec Gmsh...")
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("RectangularMeshWithAlignedLines")

        # Extraction des distances et altitudes
        distances = distance_altitude_table["Distance (m)"].to_numpy()
        z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()

        # Points du rectangle
        p1 = gmsh.model.occ.addPoint(distances[0], 0, z_coords[0])
        p2 = gmsh.model.occ.addPoint(distances[2], 0, z_coords[2])
        p3 = gmsh.model.occ.addPoint(distances[-1], 0, z_coords[-1])
        p4 = gmsh.model.occ.addPoint(distances[-2], 0, z_coords[-2])

        # Création des lignes et de la surface
        l1 = gmsh.model.occ.addLine(p1, p2)
        l2 = gmsh.model.occ.addLine(p2, p3)
        l3 = gmsh.model.occ.addLine(p3, p4)
        l4 = gmsh.model.occ.addLine(p4, p1)
        loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])
        surface = gmsh.model.occ.addPlaneSurface([loop])

        # Synchronisation du modèle
        gmsh.model.occ.synchronize()

        # Calcul des subdivisions
        length_x = max(distances) - min(distances)
        length_z = max(z_coords) - min(z_coords)
        num_div_x = max(2, math.ceil(length_x / dx_outer))
        num_div_z = max(2, math.ceil(length_z / dz_outer))

        # Ajustement des lignes verticales pour qu'elles soient alignées avec les faces des mailles
        x_line1 = 5  # Position initiale de la première ligne
        x_line2 = 6  # Position initiale de la deuxième ligne

        # Ajuster les positions pour qu'elles correspondent aux divisions des mailles
        x_line1 = distances[0] + round((x_line1 - distances[0]) / dx_outer) * dx_outer
        x_line2 = distances[0] + round((x_line2 - distances[0]) / dx_outer) * dx_outer


        # Calcul des subdivisions pour les régions
        num_div_x_outer_left = max(2, math.ceil((x_line1 - min(distances)) / dx_outer))
        num_div_x_outer_right = max(2, math.ceil((x_line2 - max(distances)) / dx_outer))
        num_div_x_inner = max(2, math.ceil((x_line2 - x_line1) / dx_inner))
        num_div_z = max(2,math.ceil(length_z / dz_outer) )



        # Création des points pour les lignes verticales
        p5 = gmsh.model.occ.addPoint(x_line1, 0,max(z_coords))
        p6 = gmsh.model.occ.addPoint(x_line1, 0, min(z_coords))
        p7 = gmsh.model.occ.addPoint(x_line2, 0, max(z_coords))
        p8 = gmsh.model.occ.addPoint(x_line2, 0, min(z_coords))

        # Création des lignes verticales
        l5 = gmsh.model.occ.addLine(p5, p6)
        l6 = gmsh.model.occ.addLine(p7, p8)

        # Application des subdivisions et recombinaison pour quadrilatères
        gmsh.model.mesh.setTransfiniteCurve(l1, num_div_x_outer_left)
        gmsh.model.mesh.setTransfiniteCurve(l3, num_div_x)
        gmsh.model.mesh.setTransfiniteCurve(l2, num_div_z)
        gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z)
        gmsh.model.mesh.setTransfiniteCurve(l5, num_div_z)
        gmsh.model.mesh.setTransfiniteCurve(l6, num_div_z)
        gmsh.model.mesh.setTransfiniteSurface(surface, "Left")
        gmsh.model.mesh.setRecombine(2, surface)

        # Génération du maillage
        gmsh.model.mesh.generate(mesh_dimension)
        gmsh.write(output_mesh_path)
        print(f"✅ Maillage sauvegardé sous : {output_mesh_path}")

    except Exception as e:
        print(f"Erreur lors de la génération du maillage : {e}")
    finally:
        gmsh.finalize()



def generate_mesh_8disc(distance_altitude_table, output_mesh_path, h_left=5,h_right=8,V_top=104.8,v_bot=103.8,  dx_outer=0.1, dz_outer=1.0, dx_inner=0.2, dz_inner=0.2,mesh_dimension=2):
    try:
        print("\nGénération du maillage structuré avec Gmsh...")
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.model.add("RectangularMesh")

        # Extraction des distances et altitudes
        distances = distance_altitude_table["Distance (m)"].to_numpy()
        z_coords = distance_altitude_table["Altitude (Z)"].to_numpy()

        # Points du rectangle
        min_z= min(z_coords)
        max_z= max(z_coords)
        p1 = gmsh.model.occ.addPoint(h_left, 0, min_z)
        p2 = gmsh.model.occ.addPoint(h_left, 0, max_z)
        p3 = gmsh.model.occ.addPoint(h_right, 0,min_z)
        p4 = gmsh.model.occ.addPoint(h_right, 0,max_z)

        # point intermediaire 
        p5 = gmsh.model.occ.addPoint(h_right, 0, v_bot)
        p6 = gmsh.model.occ.addPoint(h_left, 0, v_bot)


        # Création des lignes pour une seule surface combinée
        # Création des lignes pour une seule surface combinée
        l1 = gmsh.model.occ.addLine(p2, p4)  # Ligne supérieure
        l2 = gmsh.model.occ.addLine(p4, p5)  # Ligne droite supérieure
        l3 = gmsh.model.occ.addLine(p5, p6)  # Ligne horizontale intermédiaire (limite commune)
        l4 = gmsh.model.occ.addLine(p6, p2)  # Ligne gauche supérieure
        l5 = gmsh.model.occ.addLine(p5, p3)  # Ligne droite inférieure
        l6 = gmsh.model.occ.addLine(p3, p1)  # Ligne inférieure
        l7 = gmsh.model.occ.addLine(p1, p6)  # Ligne gauche inférieure

        # Combiner toutes les courbes dans une seule boucle
       # Création des boucles pour les deux sous-surfaces
        loop1 = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4])  # Région supérieure
        loop2 = gmsh.model.occ.addCurveLoop([l3, l5, l6, l7])  # Région inférieure

        # Création des sous-surfaces
        surface1 = gmsh.model.occ.addPlaneSurface([loop1])
        surface2 = gmsh.model.occ.addPlaneSurface([loop2])

        # Synchronisation du modèle
        gmsh.model.occ.synchronize()

        # Calcul des subdivisions
        length_x = h_right - h_left
        num_div_x = max(2, math.ceil(length_x / dx_inner))
        num_div_z_inner = max(2, math.ceil((max_z - V_top) / dz_inner))
        num_div_z_outer = max(2, math.ceil((V_top - min_z) / dz_outer))

        print(f"Nombre de divisions en X : {num_div_x}")
        print(f"Nombre de divisions en Z (intérieur) : {num_div_z_inner}")
        print(f"Nombre de divisions en Z (extérieur) : {num_div_z_outer}")

        # Application des subdivisions et recombinaison pour chaque sous-surface
        for surface in [surface1, surface2]:
            gmsh.model.mesh.setTransfiniteSurface(surface, "Left")
            gmsh.model.mesh.setRecombine(2, surface)

        # Subdivisions pour la région supérieure (surface1)
        gmsh.model.mesh.setTransfiniteCurve(l1, num_div_x)  # Ligne supérieure horizontale
        gmsh.model.mesh.setTransfiniteCurve(l3, num_div_x)  # Ligne horizontale intermédiaire (limite commune)
        gmsh.model.mesh.setTransfiniteCurve(l4, num_div_z_inner)  # Ligne gauche verticale (région supérieure)
        gmsh.model.mesh.setTransfiniteCurve(l2, num_div_z_inner)  # Ligne droite verticale (région supérieure)

        # Subdivisions pour la région inférieure (surface2)
        gmsh.model.mesh.setTransfiniteCurve(l3, num_div_x)  # Ligne horizontale intermédiaire (limite commune)
        gmsh.model.mesh.setTransfiniteCurve(l5, num_div_z_outer)  # Ligne droite verticale (région inférieure)
        gmsh.model.mesh.setTransfiniteCurve(l6, num_div_x)  # Ligne inférieure horizontale
        gmsh.model.mesh.setTransfiniteCurve(l7, num_div_z_outer)  # Ligne gauche verticale (région inférieure)

        # Génération du maillage
        gmsh.model.mesh.generate(mesh_dimension)
        gmsh.write(output_mesh_path)
        print(f"✅ Maillage sauvegardé sous : {output_mesh_path}")

    except Exception as e:
        print(f"Erreur lors de la génération du maillage : {e}")
    finally:
        gmsh.finalize()


# Appel de la fonction
output_mesh_path = os.path.join(sim_dir, "rectangular_mesh.msh")
generate_mesh_8disc(distance_altitude_table, output_mesh_path, h_left=0,h_right=16.3,V_top=104.8,v_bot=103.8,  dx_outer=0.1, dz_outer=1.0, dx_inner=0.2, dz_inner=0.05,mesh_dimension=2)

def remove_elements_above_curve(table, mesh_path):
    gmsh.initialize()
    gmsh.open(mesh_path)
    
    # Suppression de la première et dernière ligne du tableau
    table_trimmed = table.iloc[1:-1]
    d, z = table_trimmed["Distance (m)"].to_numpy(), table_trimmed["Altitude (Z)"].to_numpy()

    # Récupération des éléments du maillage (2D)
    dim = 2
    _, element_tags, node_tags = gmsh.model.mesh.getElements(dim)

    elements_to_remove = []
    for elem, nodes in zip(element_tags[0], node_tags[0].reshape(-1, 4)):
        coords = [gmsh.model.mesh.getNode(n)[0] for n in nodes]
        center_x = np.mean([c[0] for c in coords])
        center_z = np.mean([c[2] for c in coords])

        # Interpolation pour trouver l'altitude de la courbe en center_x
        z_curve = np.interp(center_x, d, z)

        # Vérification si l'élément est au-dessus de la courbe
        if center_z > z_curve:
            elements_to_remove.append(elem)

    if elements_to_remove:
        elements_to_remove = np.array(elements_to_remove, dtype=np.int32).flatten()

        # ✅ Récupérer les entités géométriques associées aux éléments (ex: la surface)
        entity_tags = gmsh.model.getEntities(dim)
        
        if entity_tags:
            tag = entity_tags[0][1]  # Prendre la première entité trouvée (surface)
            gmsh.model.mesh.removeElements(dim, tag, elements_to_remove)

            # 🔄 Reclassification des nœuds après suppression
            gmsh.model.mesh.reclassifyNodes()

            # Sauvegarde du maillage modifié
            modified_mesh_path = mesh_path.replace(".msh", "_modified.msh")
            gmsh.write(modified_mesh_path)
            print(f"✅ Modified mesh saved to: {modified_mesh_path}")
        else:
            print("⚠️ No valid entity found to remove elements from.")

    gmsh.finalize()



output_mesh_path = os.path.join(sim_dir, "rectangular_mesh.msh")
remove_elements_above_curve(distance_altitude_table, output_mesh_path)
def plot_gmsh_mesh(mesh_path):
    try:
        mesh = pv.read(mesh_path)
        plotter = pv.Plotter()
        plotter.add_mesh(mesh, show_edges=True, color="lightblue", label="Maillage Gmsh")
        points = mesh.points
# add axis coordinates
        plotter.add_axes()
#add grid
        plotter.show_grid()
        plotter.add_legend()
# add buttom zoom
        plotter.camera.zoom(1.5)



        plotter.show()
    except Exception as e:
        print(f"Erreur lors de la visualisation du maillage : {e}")

plot_gmsh_mesh(output_mesh_path.replace(".msh", "_modified.msh"))

