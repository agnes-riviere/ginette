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

# Import the original functions that are still needed
from mesh_generator import (
    plot_gmsh_mesh, process_distance_altitude, subdivisions_autour,
    creation_E_zone, remove_elements_above_curve_all_entities, readGmsh,
    voisin_mesh, maille_limite, coord_to_row_column, id_mesh_river, am_bm_modif
)

#verify if gmsh is installed
try:
    import gmsh 
except ImportError:
    print("GMSH is not installed. Please install it using 'pip install gmsh' or 'conda install gmsh'.")
    sys.exit(1)

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
            self.tag = gmsh.model.occ.addLine(p1_tag, p2_tag)
        return self.tag

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

def calculate_subdivisions(length, target_size):
    """Calculate optimal number of subdivisions for a given length and target size"""
    return max(2, math.ceil(length / target_size))

def calculate_hobo_parameters(x_hobo, x_RG, x_RD, num_div_x_centre, dx_hobo):
    """Calculate parameters for HOBO mesh refinement"""
    dx_reel = (x_RD - x_RG) / num_div_x_centre
    
    # Find which mesh cell contains the HOBO
    maille_hobo = 0
    for i in range(num_div_x_centre):
        if (x_RG + dx_reel * (i + 1)) > x_hobo:
            break
        maille_hobo += 1
    
    # Determine number of cells to add
    nb_case_a_ajouter = 2 if (x_RG + dx_reel * (maille_hobo + 1)) - x_hobo < 0.005 else 1
    
    # Calculate boundaries
    x_gauche_hobo = x_RG + dx_reel * maille_hobo
    x_droite_hobo = x_RG + dx_reel * (maille_hobo + nb_case_a_ajouter)
    
    # Calculate subdivisions
    length_x_hobo = x_droite_hobo - x_gauche_hobo
    num_div_x_hobo = calculate_subdivisions(length_x_hobo, dx_hobo)
    
    return {
        'maille_hobo': maille_hobo,
        'nb_case_a_ajouter': nb_case_a_ajouter,
        'x_gauche_hobo': x_gauche_hobo,
        'x_droite_hobo': x_droite_hobo,
        'num_div_x_hobo': num_div_x_hobo
    }

def generate_mesh_8_region_optimized(distance_altitude_table, output_mesh_path, v_bot=103.8, 
                                   x_RG=5, x_RD=17, z_riv=106, 
                                   dx_grossier=0.5, dx_precis=0.1, 
                                   hobo_points=None, dx_hobo=0.01, 
                                   dz_grossier=1.0, dz_precis=0.2, dz_hobo1=0.1, dz_hobo=0.15, 
                                   mesh_dimension=2, verbose=True,
                                   # Backward compatibility parameters
                                   x_hobo_1=None, z_hobo_1=None, x_hobo_2=None, z_hobo_2=None,verbose= True):
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
    EX : with 2 hobo

            z ↑
            |
            | p4----l103----p8----l103----p12----l103----p16
            |  |             |             |             |
            | l3           l6           l9           l12
            |  |             |             |             |
            | p3----l102----p7----l106----p11----l114----p15
            |  |             |             |             |
            | l2           l5           l8           l11
            |  |             |             |             |
            | p2----l101----p6----l105----p10----l113----p14
            |  |             |             |             |
            | l1           l4           l7           l10
            |  |             |             |             |
            | p1----l100----p5----l104----p9----l112----p13
            +-----------------------------------------------→ x
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
        z_riv = max(z_hobo_1, z_hobo_2, z_riv)
        
        if verbose:
            print(f"Domain: X=[{min_x:.2f}, {max_x:.2f}], Z=[{min_z:.2f}, {max_z:.2f}]")
            print(f"River level: {z_riv:.2f}, Bottom level: {v_bot:.2f}")

        # Calculate basic subdivisions
        length_x_gauche = x_RG - min_x
        length_x_centre = x_RD - x_RG
        length_x_droite = max_x - x_RD
        
        num_div_x_gauche = calculate_subdivisions(length_x_gauche, dx_grossier)
        num_div_x_centre = calculate_subdivisions(length_x_centre, dx_precis)
        num_div_x_droite = calculate_subdivisions(length_x_droite, dx_grossier)
        
        num_div_z_haut = calculate_subdivisions(max_z - z_riv, dz_grossier)
        num_div_z_centre = calculate_subdivisions(z_riv - v_bot, dz_precis)
        num_div_z_bas = calculate_subdivisions(v_bot - min_z, dz_grossier)

        # Calculate HOBO parameters
        if x_hobo_1 is not None and z_hobo_1 is not None:
            hobo1_params = calculate_hobo_parameters(x_hobo_1, x_RG, x_RD, num_div_x_centre, dx_hobo)
            if x_hobo_2 is not None and z_hobo_2 is not None:
                # Both HOBO points exist
                hobo2_params = calculate_hobo_parameters(x_hobo_2, x_RG, x_RD, num_div_x_centre, dx_hobo)
                num_div_x_RG_hobo1 = hobo1_params['maille_hobo']
                num_div_x_hobo1_hobo2 = calculate_subdivisions(
                    hobo2_params['x_gauche_hobo'] - hobo1_params['x_droite_hobo'], 
                    (x_RD - x_RG) / num_div_x_centre
                )
                num_div_x_hobo2_RD = num_div_x_centre - (hobo2_params['maille_hobo'] + hobo2_params['nb_case_a_ajouter'])
            else:
                # Only one HOBO point exists
                num_div_x_RG_hobo1 = hobo1_params['maille_hobo']
                num_div_x_hobo1_hobo2 = 0
                num_div_x_hobo2_RD = 0
        else:
            # No HOBO points exist
            hobo1_params = None
            hobo2_params = None
            num_div_x_RG_hobo1 = 0
            num_div_x_hobo1_hobo2 = 0
            num_div_x_hobo2_RD = 0

        if verbose:
            print(f"Subdivisions - Left: {num_div_x_gauche}, Center: {num_div_x_centre}, Right: {num_div_x_droite}")
            print(f"Z subdivisions - Top: {num_div_z_haut}, Center: {num_div_z_centre}, Bottom: {num_div_z_bas}")

        # Create all mesh points using the helper class
        points = {}

# ... (rest of the code remains the same)
        
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
        
        if hobo1_params is not None:
            points['p17'] = MeshPoint(hobo1_params['x_gauche_hobo'], 0, z_riv)
            points['p18'] = MeshPoint(hobo1_params['x_droite_hobo'], 0, z_riv)
            points['p19'] = MeshPoint(hobo1_params['x_gauche_hobo'], 0, v_bot)
            points['p20'] = MeshPoint(hobo1_params['x_droite_hobo'], 0, v_bot)
            points['p21'] = MeshPoint(hobo1_params['x_gauche_hobo'], 0, min_z)
            points['p22'] = MeshPoint(hobo1_params['x_droite_hobo'], 0, min_z)

        if hobo2_params is not None:
            points['p23'] = MeshPoint(hobo2_params['x_gauche_hobo'], 0, z_riv)
            points['p24'] = MeshPoint(hobo2_params['x_droite_hobo'], 0, z_riv)
            points['p25'] = MeshPoint(hobo2_params['x_gauche_hobo'], 0, v_bot)
            points['p26'] = MeshPoint(hobo2_params['x_droite_hobo'], 0, v_bot)
            points['p27'] = MeshPoint(hobo2_params['x_gauche_hobo'], 0, min_z)
            points['p28'] = MeshPoint(hobo2_params['x_droite_hobo'], 0, min_z)

        # Create lines using helper class
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
        if hobo1_params is not None:
            lines['l13'] = MeshLine(points['p21'], points['p19'])
            lines['l14'] = MeshLine(points['p19'], points['p17'])
            lines['l15'] = MeshLine(points['p22'], points['p20'])
            lines['l16'] = MeshLine(points['p20'], points['p18'])

        if hobo2_params is not None:
            lines['l17'] = MeshLine(points['p27'], points['p25'])
            lines['l18'] = MeshLine(points['p25'], points['p23'])
            lines['l19'] = MeshLine(points['p28'], points['p26'])
            lines['l20'] = MeshLine(points['p26'], points['p24'])

        # Horizontal lines
        lines['l100'] = MeshLine(points['p1'], points['p5'])
        lines['l101'] = MeshLine(points['p2'], points['p6'])
        lines['l102'] = MeshLine(points['p3'], points['p7'])
        lines['l103'] = MeshLine(points['p4'], points['p8'])

        # HOBO horizontal lines
        if hobo1_params is not None:
            lines['l104'] = MeshLine(points['p5'], points['p21'])
            lines['l105'] = MeshLine(points['p6'], points['p19'])
            lines['l106'] = MeshLine(points['p7'], points['p17'])
            lines['l107'] = MeshLine(points['p21'], points['p22'])
            lines['l108'] = MeshLine(points['p17'], points['p18'])
            lines['l116'] = MeshLine(points['p19'], points['p20'])

        if hobo2_params is not None:
            lines['l109'] = MeshLine(points['p22'], points['p27'])
            lines['l110'] = MeshLine(points['p20'], points['p25'])
            lines['l111'] = MeshLine(points['p18'], points['p23'])
            lines['l117'] = MeshLine(points['p27'], points['p28'])
            lines['l118'] = MeshLine(points['p25'], points['p26'])
            lines['l119'] = MeshLine(points['p23'], points['p24'])

        # Additional horizontal lines
        lines['l112'] = MeshLine(points['p9'], points['p13'])
        lines['l113'] = MeshLine(points['p10'], points['p14'])
        lines['l114'] = MeshLine(points['p11'], points['p15'])
        lines['l115'] = MeshLine(points['p12'], points['p16'])
        lines['l120'] = MeshLine(points['p28'], points['p9'])
        lines['l121'] = MeshLine(points['p26'], points['p10'])
        lines['l122'] = MeshLine(points['p24'], points['p11'])# Create lines using helper class
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
        if hobo1_params is not None:
            lines['l13'] = MeshLine(points['p21'], points['p19'])
            lines['l14'] = MeshLine(points['p19'], points['p17'])
            lines['l15'] = MeshLine(points['p22'], points['p20'])
            lines['l16'] = MeshLine(points['p20'], points['p18'])

        if hobo2_params is not None:
            lines['l17'] = MeshLine(points['p27'], points['p25'])
            lines['l18'] = MeshLine(points['p25'], points['p23'])
            lines['l19'] = MeshLine(points['p28'], points['p26'])
            lines['l20'] = MeshLine(points['p26'], points['p24'])

        # Horizontal lines
        lines['l100'] = MeshLine(points['p1'], points['p5'])
        lines['l101'] = MeshLine(points['p2'], points['p6'])
        lines['l102'] = MeshLine(points['p3'], points['p7'])
        lines['l103'] = MeshLine(points['p4'], points['p8'])

        # HOBO horizontal lines
        if hobo1_params is not None:
            lines['l104'] = MeshLine(points['p5'], points['p21'])
            lines['l105'] = MeshLine(points['p6'], points['p19'])
            lines['l106'] = MeshLine(points['p7'], points['p17'])
            lines['l107'] = MeshLine(points['p21'], points['p22'])
            lines['l108'] = MeshLine(points['p17'], points['p18'])
            lines['l116'] = MeshLine(points['p19'], points['p20'])

        if hobo2_params is not None:
            lines['l109'] = MeshLine(points['p22'], points['p27'])
            lines['l110'] = MeshLine(points['p20'], points['p25'])
            lines['l111'] = MeshLine(points['p18'], points['p23'])
            lines['l117'] = MeshLine(points['p27'], points['p28'])
            lines['l118'] = MeshLine(points['p25'], points['p26'])
            lines['l119'] = MeshLine(points['p23'], points['p24'])

        # Additional horizontal lines
        lines['l112'] = MeshLine(points['p9'], points['p13'])
        lines['l113'] = MeshLine(points['p10'], points['p14'])
        lines['l114'] = MeshLine(points['p11'], points['p15'])
        lines['l115'] = MeshLine(points['p12'], points['p16'])
        lines['l120'] = MeshLine(points['p28'], points['p9'])
        lines['l121'] = MeshLine(points['p26'], points['p10'])
        lines['l122'] = MeshLine(points['p24'], points['p11'])

        # Create regions with their respective subdivisions
        regions = []

        # Define each region with its boundary lines and subdivisions
        region_definitions = [
            ("Region1", [lines['l101'], lines['l4'], lines['l100'], lines['l1']], num_div_x_gauche, num_div_z_bas),
            ("Region2", [lines['l105'], lines['l13'], lines['l104'], lines['l4']], num_div_x_RG_hobo1, num_div_z_bas),
            ("Region3", [lines['l116'], lines['l15'], lines['l107'], lines['l13']], hobo1_params['num_div_x_hobo'], num_div_z_bas),
            ("Region4", [lines['l110'], lines['l17'], lines['l109'], lines['l15']], num_div_x_hobo1_hobo2, num_div_z_bas),
            ("Region5", [lines['l113'], lines['l10'], lines['l112'], lines['l7']], num_div_x_droite, num_div_z_bas),
            ("Region6", [lines['l102'], lines['l5'], lines['l101'], lines['l2']], num_div_x_gauche, num_div_z_centre),
            ("Region7", [lines['l106'], lines['l14'], lines['l105'], lines['l5']], num_div_x_RG_hobo1, num_div_z_centre),
            ("Region8", [lines['l111'], lines['l18'], lines['l110'], lines['l16']], num_div_x_hobo1_hobo2, num_div_z_centre),
            ("Region9", [lines['l114'], lines['l11'], lines['l113'], lines['l8']], num_div_x_droite, num_div_z_centre),
            ("Region10", [lines['l103'], lines['l6'], lines['l102'], lines['l3']], num_div_x_gauche, num_div_z_haut),
            ("Region11", [lines['l115'], lines['l12'], lines['l114'], lines['l9']], num_div_x_droite, num_div_z_haut),
        ]

        if hobo1_params is not None:
            region_definitions.append(("Region12", [lines['l108'], lines['l16'], lines['l116'], lines['l14']], hobo1_params['num_div_x_hobo'], num_div_z_centre))

        if hobo2_params is not None:
            region_definitions.append(("Region13", [lines['l118'], lines['l19'], lines['l117'], lines['l17']], hobo2_params['num_div_x_hobo'], num_div_z_bas))
            region_definitions.append(("Region14", [lines['l121'], lines['l7'], lines['l120'], lines['l19']], num_div_x_hobo2_RD, num_div_z_bas))
            region_definitions.append(("Region15", [lines['l119'], lines['l20'], lines['l118'], lines['l18']], hobo2_params['num_div_x_hobo'], num_div_z_centre))
            region_definitions.append(("Region16", [lines['l122'], lines['l8'], lines['l121'], lines['l20']], num_div_x_hobo2_RD, num_div_z_centre))

        
        # Create regions
        for region_definition in region_definitions:
            region_name, boundary_lines, num_div_x, num_div_z = region_definition
            regions.append(MeshRegion(region_name, boundary_lines, num_div_x, num_div_z))

    # Generate the mesh
        
        # Create all regions
        for name, region_lines, sub_x, sub_z in region_definitions:
            region = MeshRegion(name, region_lines, sub_x, sub_z)
            regions.append(region)

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
        gmsh.write(output_mesh_path)
        
        if verbose:
            print(f" Optimized esh saved to: {output_mesh_path}")
            print(f"Total regions created: {len(regions)}")

    except Exception as e:
        print(f" Error during mesh generation: {e}")
        raise
    finally:
        gmsh.finalize()

# Backward compatibility function
def generate_mesh_8_region(distance_altitude_table, output_mesh_path, **kwargs):
    """
    Backward compatibility wrapper for the optimized function.
    """
    return generate_mesh_8_region_optimized(distance_altitude_table, output_mesh_path, **kwargs)
