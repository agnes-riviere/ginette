#!/usr/bin/env python3
"""
Script Python pour comparer les résultats de gel et dégel 
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

def main():
    # Définition des chemins
    # path_anal = "~/run/media/ariviere/Commun/INTERFROST/Neuman-fichiers-anne/script_R_neuman/"
    # Répertoire de travail
    script_dir = os.path.expanduser("~/Programmes/ginette/application/Interfrost/Neuman/script_R_neuman")
    os.chdir(script_dir)
    
    path_gel = "../test-neuman-gel/"
    path_degel = "../test-neuman/"
    
    # Paramètres du graphique
    xlab = "temps"
    ylab = "Profondeur [m]"
    xlim = [0, 3000]
    ylim = [0, 0.1]
    
    # Lecture des données CSV
    try:
        D_anal = pd.read_csv('sol_anal.csv', header=None, sep=';', decimal='.')
        D_gel = pd.read_csv("/home/ariviere/Programmes/ginette/application/Interfrost/Neuman/test-neuman-gel/S_bound_permaf_1_t.csv", header=None, sep=';', decimal='.')
        D_degel = pd.read_csv("/home/ariviere/Programmes/ginette/application/Interfrost/Interfrost_Neuman/Conduction/S_bound_permaf_1_t.dat", header=None, sep='\s+', decimal='.')
    except FileNotFoundError as e:
        print(f"Erreur: Fichier non trouvé - {e}")
        return
    except Exception as e:
        print(f"Erreur lors de la lecture des fichiers: {e}")
        return
    
    # Configuration du graphique
    plt.figure(figsize=(5, 5), dpi=150)
    plt.margins(0)  # équivalent de xaxs="i", yaxs="i"
    
    # Tracé des courbes
    try:
        # Ligne analytique (orange)
        plt.plot(D_anal.iloc[:, 1] , D_anal.iloc[:, 3], 
                 color='orange', linewidth=1, label='Analytique')
        
        # Ligne gel (bleue)
        plt.plot(D_gel.iloc[:, 0], 2 - D_gel.iloc[:, 3], 
                 color='blue', linewidth=1, label='previous')
        
        # Ligne dégel (rouge)
        plt.plot(D_degel.iloc[:, 0], 2 - D_degel.iloc[:, 3], 
                 color='red', linewidth=1, label='ginette_V2')
        
    except Exception as e:
        print(f"Erreur lors du tracé: {e}")
        return
    
    # Configuration des axes et labels
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(xlab, fontsize=10)
    plt.ylabel(ylab, fontsize=10)
    
    # Ajustement des marges (équivalent de par(mar=c(2,2,0.5,0)))
    plt.tight_layout()
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    
    # Légende
    plt.legend(loc='upper right', fontsize=8, frameon=False)
    
    # Sauvegarde du graphique
    plt.savefig("Compar_gel_degel.png", dpi=600, bbox_inches='tight')
    
    # Affichage
    plt.show()
    
    print("Graphique généré avec succès: Compar_gel_degel.png")

if __name__ == "__main__":
    main()
