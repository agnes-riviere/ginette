#!/usr/bin/env python3
"""
Calculate volumetric calorific capacity from thermal parameters and grid search results
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def read_thermal_parameters(file_path):
    """
    Read thermal parameters from E_p_therm.dat file
    Returns dictionary with thermal properties
    """
    thermal_params = {}
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        if 'cpe=' in line:
            # Specific heat capacity of water (J/kg/K)
            thermal_params['cpe'] = float(line.split('=')[1].split()[0].replace('D', 'E'))
        elif 'cpm=' in line:
            # Specific heat capacity of solid (J/kg/K)  
            thermal_params['cpm'] = float(line.split('=')[1].split()[0])
        elif 'alandam=' in line:
            # Thermal conductivity of solid (W/m/K)
            thermal_params['alandam'] = float(line.split('=')[1].split()[0].replace('D', 'E'))
        elif 'alandae=' in line:
            # Thermal conductivity of water (W/m/K)
            thermal_params['alandae'] = float(line.split('=')[1].split()[0].replace('D', 'E'))
    
    return thermal_params

def read_results_data(file_path):
    """
    Read grid search results
    """
    return pd.read_csv(file_path, sep=' ')

def calculate_volumetric_calorific_capacity(thermal_params, porosity, rhos):
    """
    Calculate volumetric calorific capacity for one line
    
    Parameters:
    - thermal_params: dict with thermal parameters
    - porosity: porosity value for this line
    - rhos: solid density value for this line
    
    Returns:
    - volumetric calorific capacity (J/m³/K) for this line
    """
    # Extract parameters from thermal parameter file (must be present)
    if 'cpe' not in thermal_params:
        raise ValueError("cpe (water specific heat capacity) not found in thermal parameters file")
    if 'cpm' not in thermal_params:
        raise ValueError("cpm (solid specific heat capacity) not found in thermal parameters file")
    
    cpe = thermal_params['cpe']    # J/kg/K - water specific heat capacity
    cpm = thermal_params['cpm']    # J/kg/K - solid specific heat capacity
    
    rho_water = 1000  # kg/m³ - water density
    
    # Volumetric calorific capacity calculation for this line
    # C_vol = φ * ρ_water * c_water + (1-φ) * ρ_solid * c_solid
    C_vol = porosity * rho_water * cpe + (1 - porosity) * rhos * cpm
    
    return C_vol

def main():
    """
    Main function to calculate volumetric calorific capacity
    """
    print("=" * 60)
    print("VOLUMETRIC CALORIFIC CAPACITY CALCULATION")
    print("=" * 60)
    
    # File paths
    thermal_file = './SYNTHETIC_CASES/E_p_therm.dat'
    results_file = './results/results.txt'
    
    try:
        # Read thermal parameters
        print("\n1. Reading thermal parameters...")
        thermal_params = read_thermal_parameters(thermal_file)
        print(f"   Thermal parameters loaded:")
        for key, value in thermal_params.items():
            print(f"   - {key}: {value}")
        
        # Read results data
        print("\n2. Reading results data...")
        results_df = read_results_data(results_file)
        
        # Convert columns to numeric types to avoid string/sequence errors
        results_df['poro'] = pd.to_numeric(results_df['poro'], errors='coerce')
        results_df['cap'] = pd.to_numeric(results_df['cap'], errors='coerce')
        
        print(f"   Grid search results loaded: {len(results_df)} simulations")
        print(f"   Porosity range: {results_df['poro'].min():.3f} to {results_df['poro'].max():.3f}")
        print(f"   Cap range: {results_df['cap'].min():.1f} to {results_df['cap'].max():.1f}")
        
        # Calculate volumetric calorific capacity
        print("\n3. Calculating volumetric calorific capacity...")
        print(f"   Read cpe (water specific heat): {thermal_params['cpe']} J/kg/K")
        print(f"   Read cpm (solid specific heat): {thermal_params['cpm']} J/kg/K") 
        print(f"   Using rho_water: 1000 kg/m³")
        
        # Create results summary - calculate for each line
        results_summary = []
        
        for i in range(len(results_df)):
            # Get values for this line and ensure they are numeric
            porosity = float(results_df.iloc[i]['poro'])
            rhos = float(results_df.iloc[i]['cap'])
            
            # Calculate C_vol for this line
            C_vol = calculate_volumetric_calorific_capacity(thermal_params, porosity, rhos)
            
            results_summary.append({
                'ID': results_df.iloc[i]['ID'],
                'porosity': porosity,
                'cap': rhos,
                'C_volumetric_J_m3_K': C_vol,
                'C_volumetric_MJ_m3_K': C_vol / 1e6  # Convert to MJ/m³/K
            })
            
            # Print calculation details for first 5 lines
            if i < 5:
                print(f"   Ligne {i+1}: Porosity={porosity:.3f}, Cap={rhos:.1f} -> C_vol={C_vol:.0f} J/m³/K ({C_vol/1e6:.3f} MJ/m³/K)")
        
        print(f"   ... Calcul terminé pour {len(results_df)} lignes")
        
        # Create summary DataFrame
        summary_df = pd.DataFrame(results_summary)
        
        # Save results
        output_file = './results/volumetric_calorific_capacity.csv'
        summary_df.to_csv(output_file, index=False)
        print(f"\n4. Results saved to: {output_file}")
        
        # Plot results
        plt.figure(figsize=(12, 8))
        
        # Subplot 1: C_vol vs Porosity
        plt.subplot(2, 1, 1)
        plt.scatter(summary_df['porosity'], summary_df['C_volumetric_MJ_m3_K'], 
                   c=summary_df['cap'], cmap='viridis', alpha=0.7)
        plt.colorbar(label='Cap')
        plt.xlabel('Porosity')
        plt.ylabel('Volumetric Calorific Capacity (MJ/m³/K)')
        plt.title('Volumetric Calorific Capacity vs Porosity (colored by Cap)')
        plt.grid(True, alpha=0.3)
        
        # Subplot 2: C_vol vs Cap
        plt.subplot(2, 1, 2)
        plt.scatter(summary_df['cap'], summary_df['C_volumetric_MJ_m3_K'], 
                   c=summary_df['porosity'], cmap='plasma', alpha=0.7)
        plt.colorbar(label='Porosity')
        plt.xlabel('Cap')
        plt.ylabel('Volumetric Calorific Capacity (MJ/m³/K)')
        plt.title('Volumetric Calorific Capacity vs Cap (colored by Porosity)')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        plot_file = './results/volumetric_calorific_capacity.png'
        plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        print(f"5. Plot saved to: {plot_file}")
        plt.show()
        
        print("\n" + "=" * 60)
        print("CALCULATION COMPLETED SUCCESSFULLY")
        print("=" * 60)
        
        return summary_df
        
    except FileNotFoundError as e:
        print(f"\nFile not found: {e}")
        print("Please check file paths.")
    except Exception as e:
        print(f"\nError occurred: {e}")
        print("Please check file paths and data format.")

if __name__ == "__main__":
    main()

