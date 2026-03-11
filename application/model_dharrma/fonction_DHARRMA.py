import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.cm import copper

dossier_actuel = Path(__file__).parent

def smooth_square_wave(x, period):
    return 2 / (1 + np.exp(-5 * np.sin(2 * np.pi * x / period))) - 1

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def creation_temp(jours_simulation,pas_simulation,bottom_temp,deg_per_day,dayli_fluctuation,weekly_fluctuation,temp_offset):
    '''
    
    '''
    df_BC = pd.DataFrame()
    df_BC['times'] = np.arange(0, jours_simulation*86400+pas_simulation, pas_simulation)
    df_BC['days'] = df_BC['times'] /  (60*60*24)
    # Calculate the top boundary temperature
    df_BC['T_top'] = (dayli_fluctuation * np.sin(2 * np.pi * df_BC['times'] / (60*60*24))  # Daily fluctuation
                  + deg_per_day * df_BC['days']  # Linear increase per day
                  + weekly_fluctuation * np.sin(2 * np.pi * df_BC['times'] / (60*60*24*7))  # Weekly fluctuation
                  + temp_offset)  # Base temperature offset
    
    #Set the bottom boundary temperature to a constant value 
    df_BC['T_bottom'] = bottom_temp * np.ones_like(df_BC['times'])

    df_BC[['T_top','T_bottom']].to_csv(f"{dossier_actuel}/input_ginette/E_temp_t.dat", sep=" ", index = False,header=False)

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def creation_infiltration(jours_simulation,pas_simulation,min_infiltration,max_infiltration,duree_transition,debut_jour_infiltration,duree_max_infiltration):
    # nbr_tot = 86400 * jours_simulation/pas_simulation
    infiltration = list()
    cumul_list = list()
    cumul = 0

    coeff_directeur_transition = (max_infiltration-min_infiltration)/duree_transition

    i_infiltration = 0

    for i in range (900,jours_simulation*86400,pas_simulation):
    # temps_sec = i*60

        temps_jour = i/86400
        if i_infiltration == len(debut_jour_infiltration):
            valeur_infiltration = min_infiltration
            # infiltration.append(0)
            # time.append(temps_jour)

        elif temps_jour <= debut_jour_infiltration[i_infiltration]:
            valeur_infiltration = min_infiltration
            # infiltration.append(0)
            # time.append(temps_jour)
        
        elif temps_jour >= debut_jour_infiltration[i_infiltration]:
            if temps_jour >= debut_jour_infiltration[i_infiltration] + duree_transition*2 + duree_max_infiltration : # Fin de scénarion
                i_infiltration = i_infiltration +1
                valeur_infiltration = min_infiltration
                # infiltration.append(0)
                # time.append(temps_jour)

            elif temps_jour > debut_jour_infiltration[i_infiltration] + duree_transition + duree_max_infiltration : #Pente descendante
                valeur_infiltration = -coeff_directeur_transition*(temps_jour-debut_jour_infiltration[i_infiltration]-duree_transition-duree_max_infiltration) + max_infiltration
                # infiltration.append(valeur_infiltration)
                # time.append(temps_jour)

            elif temps_jour > debut_jour_infiltration[i_infiltration] + duree_transition : # Dans le max
                valeur_infiltration = max_infiltration
                # infiltration.append(max_infiltration)
                # time.append(temps_jour)

            else : # Pente ascendante
                valeur_infiltration = coeff_directeur_transition*(temps_jour-debut_jour_infiltration[i_infiltration])
        
        cumul = cumul + valeur_infiltration
        cumul_list.append(cumul)
        infiltration.append(valeur_infiltration)
        # time.append(temps_jour)

    # Ecriture pluie
    with open(f"{dossier_actuel}/input_ginette/E_debit_haut_t.dat", "w") as f:
        for item in infiltration:
            f.write(str(item) + "\n") 

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def creation_infiltration_evapo(jours_simulation,pas_simulation,min_infiltration,max_infiltration,duree_transition_infiltration,debut_jour_infiltration,duree_max_infiltration,min_evapo,max_evapo,duree_transition_evapo,debut_jour_evapo,duree_max_evapo):
    # nbr_tot = 86400 * jours_simulation/pas_simulation
    infiltration = list()

    # ### Test chevauchement infiltration et evapo
    # for debut_infiltration in debut_jour_infiltration :
    #     fin_infiltration = debut_infiltration + duree_max_infiltration + 2 * duree_transition_infiltration
    #     for debut_evapo in debut_jour_evapo :
    #         fin_evapo_transpiration = debut_evapo + duree_max_evapo + 2 * duree_transition_evapo

    #         if 


    coeff_directeur_transition = (max_infiltration-min_infiltration)/duree_transition_infiltration
    coeff_directeur_transition_evapo = -(max_evapo-min_evapo)/duree_transition_evapo

    i_infiltration = 0
    i_evapo = 0

    for i in range (900,jours_simulation*86400,pas_simulation):
    # temps_sec = i*60

        temps_jour = i/86400
        if i_infiltration == len(debut_jour_infiltration):
            
            if i_evapo != len(debut_jour_evapo) :
                if temps_jour <= debut_jour_evapo[i_evapo]:
                    valeur_infiltration = min_infiltration
                elif temps_jour >= debut_jour_evapo[i_evapo] + duree_transition_evapo*2 + duree_max_evapo : # Fin de scénarion
                    i_evapo = i_evapo +1
                    valeur_infiltration = min_infiltration

                elif temps_jour > debut_jour_evapo[i_evapo] + duree_transition_evapo + duree_max_evapo : #Pente descendante
                    valeur_infiltration = -coeff_directeur_transition_evapo*(temps_jour-debut_jour_evapo[i_evapo]-duree_transition_evapo-duree_max_evapo) - max_evapo

                elif temps_jour > debut_jour_evapo[i_evapo] + duree_transition_evapo : # Dans le max
                    valeur_infiltration = -max_evapo


                else : # Pente ascendante
                    valeur_infiltration = coeff_directeur_transition_evapo*(temps_jour-debut_jour_evapo[i_evapo])
            
            else :
                valeur_infiltration = min_infiltration

        elif temps_jour <= debut_jour_infiltration[i_infiltration]:
            # valeur_infiltration = min_infiltration
            # AJOUTER EVAPO Transpiration ici
            # print(i_evapo)
            if i_evapo != len(debut_jour_evapo):
                if temps_jour <= debut_jour_evapo[i_evapo]:
                    valeur_infiltration = min_infiltration
                elif temps_jour >= debut_jour_evapo[i_evapo] + duree_transition_evapo*2 + duree_max_evapo : # Fin de scénarion
                    i_evapo = i_evapo +1
                    valeur_infiltration = min_infiltration

                elif temps_jour > debut_jour_evapo[i_evapo] + duree_transition_evapo + duree_max_evapo : #Pente descendante
                    valeur_infiltration = -coeff_directeur_transition_evapo*(temps_jour-debut_jour_evapo[i_evapo]-duree_transition_evapo-duree_max_evapo) - max_evapo

                elif temps_jour > debut_jour_evapo[i_evapo] + duree_transition_evapo : # Dans le max
                    valeur_infiltration = -max_evapo


                else : # Pente ascendante
                    valeur_infiltration = coeff_directeur_transition_evapo*(temps_jour-debut_jour_evapo[i_evapo])
        
            else :
                valeur_infiltration = min_infiltration

        elif temps_jour >= debut_jour_infiltration[i_infiltration]:
            if temps_jour >= debut_jour_infiltration[i_infiltration] + duree_transition_infiltration*2 + duree_max_infiltration : # Fin de scénarion
                i_infiltration = i_infiltration +1
                valeur_infiltration = min_infiltration
                # infiltration.append(0)
                # time.append(temps_jour)

            elif temps_jour > debut_jour_infiltration[i_infiltration] + duree_transition_infiltration + duree_max_infiltration : #Pente descendante
                valeur_infiltration = -coeff_directeur_transition*(temps_jour-debut_jour_infiltration[i_infiltration]-duree_transition_infiltration-duree_max_infiltration) + max_infiltration
                # infiltration.append(valeur_infiltration)
                # time.append(temps_jour)

            elif temps_jour > debut_jour_infiltration[i_infiltration] + duree_transition_infiltration : # Dans le max
                valeur_infiltration = max_infiltration
                # infiltration.append(max_infiltration)
                # time.append(temps_jour)

            else : # Pente ascendante
                valeur_infiltration = coeff_directeur_transition*(temps_jour-debut_jour_infiltration[i_infiltration])
        
        infiltration.append(valeur_infiltration)
        # time.append(temps_jour)

    # Ecriture pluie
    with open(f"{dossier_actuel}/input_ginette/E_debit_haut_t.dat", "w") as f:
        for item in infiltration:
            f.write(str(item) + "\n") 

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def decomposition_nombre(x):
    if '.' in str(x):
        entier = str(x).replace('.','')
        puissance = -len(str(x).split('.')[-1])
    elif ',' in str(x):
        entier = str(x).replace(',','')
        puissance = -len(str(x).split(',')[-1])
    else :
        return x,0
    return int(entier),puissance

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_graphe_pluie(path_pluie_A,jours_sim, pas_sim,path_pluie_B = None, liste_mesure = [], hauteur_WT_A = None, hauteur_WT_B = None, barre_vertical = [],path_fig =None, max_cumul = 1E-04):
    '''
    Fonction permettant de plotter le ou les scénario de pluie pour l'article avec la hauter d la WT correspondante

    Entrée :
    path_pluie_A (str) : chemin absolue ammenant à la première pluie que l'on veut plotter. Si path_pluie_B n'est pas donné seule cette pluie sera plotté.
    jours_sim (int) : Nombre de jour de simulation que l'on veut plotter.
    pas_sim (int) : Temps correspondant au pas de la simulation ATTENTION IL FAUT METTRE CE TEMPS EN SEC !!!!!!!!!!!!!!!
    path_pluie_B (str) : idem que pluie A mais lui n'est pas obligatoire permet de comparer 2 scénarios de pluie.
    list_mesure (list) : Liste permettant de mettre des barres verticales sur les graphs correspondant au mesure de profil.
    hauteur_WT_A (str) : chemin amenant au fichier du la hauteur de la WT obtenue avec la fonction depth_WT_time. Plotter si donné.
    hauteru_WT_B (str) : idem que pour la A. Plotter si donné et path_pluie_B aussi.
    barre_vertical (list) : idem que liste mesure avec une représentation différente.
    path_fig(str) : Chemin permettant de sauvgarder la figure si l'on veut    
    '''
    
    pluie_A = np.loadtxt(path_pluie_A)

    cumul_A_list = list()
    cumul_A = 0

    time = list()
    i=0

    if path_pluie_B:
        pluie_B = np.loadtxt(path_pluie_B)

        cumul_B_list = list()
        cumul_B = 0

        for temps_sec in range(900,jours_sim*86400,pas_sim):
            temps_jour = temps_sec/86400

            cumul_A = cumul_A + pluie_A[i]
            cumul_A_list.append(cumul_A)
            cumul_B = cumul_B + pluie_B[i]
            cumul_B_list.append(cumul_B)

            time.append(temps_jour)
            i=i+1

        if hauteur_WT_A and hauteur_WT_B :
            depth_WT_A = np.loadtxt(hauteur_WT_A)
            depth_WT_B = np.loadtxt(hauteur_WT_B)


            fig, ax = plt.subplots(2,2,figsize=(20,7))
            ax[0,0].plot(time,pluie_A,label = 'infiltration')
            ax[0,0].fill_between(time, pluie_A, 0, alpha=0.3, color="blue")
            ax[0,0].set_xlim(0,120)
            ax[0,0].set_ylim(0, 3.2E-08)
            ax[0,0].tick_params(axis='both', labelsize=20)
            ax[0,0].set_xticklabels([])
            # ax[0,0].set_xlabel('day',fontsize = 20)
            ax[0,0].set_ylabel('Infiltration (m/s)', fontsize = 20)

            ax2 = ax[0,0].twinx()
            ax2.plot(time,cumul_A_list,color = 'black', label = 'cumul')
            ax2.set_ylim(0, max_cumul)
            ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            ax2.tick_params(axis='both', labelsize=20)
            ax2.set_ylabel('Rain cumulation (m)', fontsize = 20)
            # Fusionner les handles (lignes) des deux axes pour une seule légende
            lines_1, labels_1 = ax[0,0].get_legend_handles_labels()
            lines_2, labels_2 = ax2.get_legend_handles_labels()
            ax[0,0].legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right")

            ax[0,1].plot(time,depth_WT_A, label = 'WT level')
            ax[0,1].tick_params(axis='both', labelsize=20)
            ax[0,1].set_ylabel('Depth (m)', fontsize = 20)
            ax[0,1].set_ylim(1.97, 2.2)
            # ax[0,1].set_xlabel("day", fontsize = 20)
            ax[0,1].set_xticklabels([])
            
            ax[1,0].plot(time,pluie_B,label = 'infiltration')
            ax[1,0].fill_between(time, pluie_B, 0, alpha=0.3, color="blue")
            ax[1,0].set_xlim(0,120)
            ax[1,0].set_ylim(0, 3.2E-08)
            ax[1,0].tick_params(axis='both', labelsize=20)
            ax[1,0].set_xlabel('day',fontsize = 20)
            ax[1,0].set_ylabel('Infiltration (m/s)', fontsize = 20)

            ax2 = ax[1,0].twinx()
            ax2.plot(time,cumul_B_list,color = 'black', label = 'cumul')
            ax2.set_ylim(0, max_cumul)
            ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            ax2.tick_params(axis='both', labelsize=20)
            ax2.set_ylabel('Rain cumulation (m)', fontsize = 20)
            lines_1, labels_1 = ax[1,0].get_legend_handles_labels()
            lines_2, labels_2 = ax2.get_legend_handles_labels()
            ax[1,0].legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right")

            ax[1,1].plot(time,depth_WT_B, label = 'WT level')
            ax[1,1].tick_params(axis='both', labelsize=20)
            ax[1,1].set_ylabel('Depth (m)', fontsize = 20)
            ax[1,1].set_ylim(1.97, 2.2)
            ax[1,1].set_xlabel("day", fontsize = 20)

        else :
            fig, ax = plt.subplots(2,1,figsize=(11,8))
            # Plot 1
            ax[0].plot(time,pluie_A,label = 'infiltration')
            ax[0].fill_between(time, pluie_A, 0, alpha=0.3, color="blue")
            ax[0].set_xlim(0,120)
            ax[0].set_ylim(0, 3.2E-08)
            ax[0].tick_params(axis='both', labelsize=20)
            # ax[0].set_xlabel('day',fontsize = 20)
            ax[0].set_xticklabels([])
            ax[0].set_ylabel('Infiltration (m/s)', fontsize = 20)

            ax2 = ax[0].twinx()
            ax2.plot(time,cumul_A_list,color = 'black', label = 'cumul')
            ax2.set_ylim(0, max_cumul)
            ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            ax2.tick_params(axis='both', labelsize=20)
            ax2.set_ylabel('Rain cumulation (m)', fontsize = 20)

            # Fusionner les handles (lignes) des deux axes pour une seule légende
            lines_1, labels_1 = ax[0].get_legend_handles_labels()
            lines_2, labels_2 = ax2.get_legend_handles_labels()
            ax[0].legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right")

            ax[1].plot(time,pluie_B,label = 'infiltration')
            ax[1].fill_between(time, pluie_B, 0, alpha=0.3, color="blue")
            ax[1].set_xlim(0,120)
            ax[1].set_ylim(0, 3.2E-08)
            ax[1].tick_params(axis='both', labelsize=20)
            ax[1].set_xlabel('day',fontsize = 20)
            ax[1].set_ylabel('Infiltration (m/s)', fontsize = 20)

            ax2 = ax[1].twinx()
            ax2.plot(time,cumul_B_list,color = 'black', label = 'cumul')
            ax2.set_ylim(0, max_cumul)
            ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            ax2.tick_params(axis='both', labelsize=20)
            ax2.set_ylabel('Rain cumulation (m)', fontsize = 20)
            lines_1, labels_1 = ax[1].get_legend_handles_labels()
            lines_2, labels_2 = ax2.get_legend_handles_labels()
            ax[1].legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right")
    

    
    else :
        for temps_sec in range(900,jours_sim*86400,pas_sim):
            temps_jour = temps_sec/86400

            cumul_A = cumul_A + pluie_A[i]
            cumul_A_list.append(cumul_A)

            time.append(temps_jour)
            i=i+1

            # print(i)
        
        if hauteur_WT_A :
            depth_WT_A = np.loadtxt(hauteur_WT_A)

            fig, ax = plt.subplots(1,2,figsize=(15,4))

            ax[0].plot(time,pluie_A,label = 'infiltration')
            ax[0].fill_between(time, pluie_A, 0, alpha=0.3, color="blue")
            ax[0].set_xlim(0,120)
            # ax[0].set_ylim(0, 3.2E-08)
            ax[0].tick_params(axis='both', labelsize=20)
            ax[0].set_xlabel('day',fontsize = 20)
            ax[0].set_ylabel('Infiltration (m/s)', fontsize = 20)

            ax2 = ax[0].twinx()
            ax2.plot(time,cumul_A_list,color = 'black', label = 'cumul')
            ax2.set_ylim(0, max_cumul)
            ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
            ax2.tick_params(axis='both', labelsize=20)
            ax2.set_ylabel('Rain cumulation (m)', fontsize = 18)

            # Fusionner les handles (lignes) des deux axes pour une seule légende
            lines_1, labels_1 = ax[0].get_legend_handles_labels()
            lines_2, labels_2 = ax2.get_legend_handles_labels()
            ax[0].legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right")

            ax[1].plot(time,depth_WT_A, label = 'WT level')
            ax[1].tick_params(axis='both', labelsize=20)
            ax[1].set_ylabel('Depth (m)', fontsize = 18)
            ax[1].set_xlabel("day", fontsize = 20)
            
        
        else :
            fig, ax = plt.subplots(figsize=(11,4))
            ax.plot(time,pluie_A,label = 'infiltration')
            ax.fill_between(time, pluie_A, 0, alpha=0.3, color="blue")
            ax2 = ax.twinx()
            ax2.plot(time,cumul_A_list,color = 'black', label = 'cumul')
            ax.set_xlim(0,120)
            # ax.set_ylim(0, 3.2E-08)
            ax2.set_ylim(0, max_cumul)
            ax.tick_params(axis='both', labelsize=20)
            ax.set_xlabel('day',fontsize = 20)
            ax.set_ylabel('Infiltration (m/s)', fontsize = 20)
            ax2.tick_params(axis='both', labelsize=20)
            ax2.set_ylabel('Rain cumulation (m)', fontsize = 20)

            # Fusionner les handles (lignes) des deux axes pour une seule légende
            lines_1, labels_1 = ax.get_legend_handles_labels()
            lines_2, labels_2 = ax2.get_legend_handles_labels()
            ax.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper right")
    
    if len(barre_vertical) != 0 :
        if path_pluie_B == None and hauteur_WT_A == None and hauteur_WT_B == None:
            for i in liste_mesure:
                ax.axvline(x=i, color="r", linestyle="--")
        else :
            for a in ax.flat:
                for i in barre_vertical:
                    a.axvline(x=i, color="r", linestyle="--")

    if len(liste_mesure) != 0 :
        if path_pluie_B == None and hauteur_WT_A == None and hauteur_WT_B == None:
            for i in liste_mesure:
                ax.axvline(x=i, color="green", linestyle="-")
        else :
            for a in ax.flat:
                for i in liste_mesure:
                    a.axvline(x=i, color="green", linestyle="-")
    
    plt.tight_layout()
    if path_fig:
        plt.savefig(path_fig)
    # plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def lorentzian_error(v_picked, f_picked, dx, Nx, a=0.5):
    # Factor to adapt error depending on window size
    fac = 10 ** (1 / np.sqrt(Nx * dx))

    # Resolution
    Dc_left = 1 / (1 / v_picked - 1 / (2 * f_picked * Nx * fac * dx))
    Dc_right = 1 / (1 / v_picked + 1 / (2 * f_picked * Nx * fac * dx))
    Dc = np.abs(Dc_left - Dc_right)

    # Absolute uncertainty
    dc = (10**-a) * Dc

    for i, (err, v) in enumerate(zip(dc, v_picked)):
        if err > 0.4 * v:
            dc[i] = 0.4 * v
        if err < 5:
            dc[i] = 5

    return dc

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def depth_WT_time (path_sat, jours_sim, pas_sim,path_hauteur):

    '''
    Fonction qui calcule la hauteur de la WT en fonction du temps puis la sauvegarde dans un fichier.
    A modifer pour le sable pour éviter de prendre en compte le front de saturation + faire attention au dernier élément
    Entrée :
    path_sat (str) : chemin absolue ammenant à la saturation par maille par pas de temps qui va permettre de calculer la hauteur de la WT
    jours_sim (int) : Nombre de jour de simulation que l'on veut plotter.
    pas_sim (int) : Temps correspondant au pas de la simulation ATTENTION IL FAUT METTRE CE TEMPS EN SEC !!!!!!!!!!!!!!!
    path_hauteur (str) : Chemin permettant de sauvgarder le fichier dde sortie comtenant la hauteur de la WT en fonction du temps.    
    '''

    saturation = np.loadtxt(path_sat)
    WT_depth = []

    avancement_mem = None

    dt, z, Sw = saturation[:,0], saturation[:,1], saturation[:,2]
    df_sat = pd.DataFrame({'dt': dt,'z': z, 'Sw': Sw})

    for t_profil in range(900,jours_sim*86400,pas_sim):

        
        df_sat_profil = df_sat[(df_sat["dt"] == t_profil)].copy()

        # Boucle pour savoir le z_sat cad la hauteur max où la saturation est égale à 1
        z_sat = None
        i=0
        while (i < len(df_sat_profil['z']) and z_sat == None) :
                
                if df_sat_profil.iloc[i]['Sw'] == 1.000 :
                    z_sat = df_sat_profil.iloc[i]['z']
                else :
                    i = i+1
        if (z_sat == None) : print ('Non saturée')

        WT_depth.append(z_sat)
        avancement = int(t_profil/10368000 *100)
        if avancement_mem != avancement :
            print(f'{avancement} %')
            avancement_mem = avancement
    with open(path_hauteur, "w") as f:
        for item in WT_depth:
            f.write(str(item) + "\n") 

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def three_plot_output_ginette_dharrma(path_result,debut_rp,jours_rp,pas_rp,path_save_fig=None, barre_vertical = [],lim_depth=2.1):
    '''
    Fonction permettant de plotter et sauvegarder les figures représentant les résultats de sortie de ginette pour l'article.


    Entrée :
    path_result (str) : chemin absolue ammenant au dossier contenant les fichiers résultats à plotter.
    pluie (str) : Lettre permettant de savoir quelle scénario de pluie on veut plotter. Inutile quand diff = True.
    debut_rp (int) : Nombre (en jour) qui correspond au début de la représentation graphique
    jours_rp (int) : Nombre de jour de simulation que l'on veut plotter.
    pas_rp (int) : Temps (en jour) correspondant au pas de représentation
    list_mesure (list) : Liste permettant de mettre des barres verticales sur les graphs correspondant au mesure de profil.
    diff (bool) : Si diff = True la fonction plot la différence entre le scénario A et B. Sinon plot uniquement la pluie séléctionné avec "pluie"
    barre_vertical (list) : idem que liste mesure avec une représentation différente.
    path_fig(str) : Chemin permettant de sauvgarder la figure si l'on veut    
    '''
    


    vit_eau = np.loadtxt(f'{path_result}/input_ginette/S_vitesse_profil.dat') # Vitesse du flux d'eau 
    pression_txt = np.loadtxt(f'{path_result}/input_ginette/S_pressure_profil_t.dat')
    sat_txt = np.loadtxt(f'{path_result}/input_ginette/S_saturation_profil_t.dat')
        
    dt = vit_eau[:,0]
    z = vit_eau[:,1]
    vit = vit_eau[:,2]
    pression = pression_txt[:,2]
    sat = sat_txt[:,2]

    temps_voulu = []

    df_tot = pd.DataFrame({'dt': dt,'z': z, 'vit': -vit, 'pression': pression, 'sat' : sat})
    # print(df_tot)
    for t_profil in range (debut_rp*86400,jours_rp*86400,pas_rp*86400):
        temps_voulu.append(t_profil)
    
    # print(temps_voulu)

    df_profil_filtre = df_tot[df_tot["z"]>=-lim_depth].copy()
    # print(df_profil_filtre)
    df_profil_filtre = df_profil_filtre[df_profil_filtre["dt"].isin(temps_voulu)].copy()
    df_profil_filtre["z"] = df_profil_filtre["z"]

        

    # Utilisation de pivot ----------------------------------------------------------------------------------------------------------------
    grid_vit = df_profil_filtre.pivot(index='z',columns='dt',values='vit')
    # print(grid_vit)
    grid_pression = df_profil_filtre.pivot(index='z',columns='dt',values='pression')
    grid_sat = df_profil_filtre.pivot(index='z',columns='dt',values='sat')
    # print("Grid shape:", grid.shape)

    # Plot --------------------------------------------------------------------------------------------------------------------------------
    plt.rcParams['font.size']=20
    fig, ax = plt.subplots(1,3, figsize=(20,5))

    im0 = ax[0].imshow(grid_vit, cmap="viridis",aspect='auto', origin='lower',extent=[0, jours_rp,
            grid_vit.index.min(), grid_vit.index.max()])
    # ax[0].colorbar()
    ax[0].set_ylabel('Depth (m)')
    ax[0].set_xlabel('Day')
    ax[0].set_title('Water velocity (m/s)')

    im1 = ax[1].imshow(grid_pression, cmap="viridis", aspect='auto', origin='lower',extent=[0, jours_rp,
            grid_vit.index.min(), grid_vit.index.max()])
    ax[1].set_yticklabels([])
    ax[1].set_xlabel('Day')
    ax[1].set_title('Presure (Pa)')

    im2 = ax[2].imshow(grid_sat, cmap="viridis", aspect='auto', origin='lower',extent=[0, jours_rp,
            grid_vit.index.min(), grid_vit.index.max()])
    ax[2].set_yticklabels([])
    ax[2].set_xlabel('Day')
    ax[2].set_title('Saturation (-)')

    fig.colorbar(im0, ax=ax[0])#, label="Velocity (m/s)")
    cbar1 = fig.colorbar(im1, ax=ax[1])#, label="Presure (Pa)")
    cbar1.ax.ticklabel_format(style='sci', scilimits=(0,0))
    fig.colorbar(im2, ax=ax[2])#, label="Saturation (-)")
    if len(barre_vertical) != 0 :
        for a in ax.flat:
            for i in barre_vertical:
                a.axvline(x=i, color="r", linestyle="--")
    plt.tight_layout()
    if path_save_fig :
        plt.savefig(path_save_fig)
    # plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def three_plot_propriete_geophy_dharrma(path_result,debut_rp,jours_rp,pas_rp, facies, barre_vertical=[],path_save_fig =None,lim_depth=2.1,Vp_Vs = 'Vs'):
    '''
    Fonction permettant de plotter et sauvegarder les figures représentant les propriété geophysique obtenu avec le FWD model pour l'article.


    Entrée :
    path_result (str) : chemin absolue ammenant au dossier contenant les fichiers résultats à plotter.
    pluie (str) : Lettre permettant de savoir quelle scénario de pluie on veut plotter. Inutile quand diff = True.
    debut_rp (int) : Nombre (en jour) qui correspond au début de la représentation graphique
    jours_rp (int) : Nombre de jour de simulation que l'on veut plotter.
    pas_rp (int) : Temps (en jour) correspondant au pas de représentation
    list_mesure (list) : Liste permettant de mettre des barres verticales sur les graphs correspondant au mesure de profil.
    diff (bool) : Si diff = True la fonction plot la différence entre le scénario A et B. Sinon plot uniquement la pluie séléctionné avec "pluie"
    barre_vertical (list) : idem que liste mesure avec une représentation différente.
    path_save_fig(str) : Chemin permettant de sauvgarder la figure si l'on veut   
    '''


    path_temp = f'{path_result}/input_ginette/S_temperature_t.dat'
    path_Vs = f'{path_result}/sismique/{facies}/output_SL_kk4_Vp_Vs.dat'
    path_rho_vrai = f'{path_result}/elec/{facies}/rho_vrai.dat'
    
    temp_txt = np.loadtxt(path_temp)
    Vs_txt = np.loadtxt(path_Vs)
    rho_vrai_txt = np.loadtxt(path_rho_vrai)

    dt_temp = temp_txt[:,0]
    z_temp = temp_txt[:,1]
    temp = temp_txt[:,2]

    dt = Vs_txt[:,0]
    z = Vs_txt[:,1]
    if Vp_Vs == 'Vp':
        Vs = Vs_txt[:,2]
    else :
        Vs = Vs_txt[:,3]
    rho_vrai = rho_vrai_txt[:,2]

    temps_voulu = []

    df_tot_Vs_rho = pd.DataFrame({'dt': dt,'z': z, 'Vs': Vs, 'rho_vrai' : rho_vrai})
    df_tot_temp = pd.DataFrame({'dt_temp': dt_temp,'z_temp':z_temp, 'temp': temp})

    for t_profil in range (debut_rp*86400,jours_rp*86400,pas_rp*86400):
            temps_voulu.append(t_profil)

    df_profil_filtre_Vs_rho = df_tot_Vs_rho[df_tot_Vs_rho["z"]>=-lim_depth].copy()#
    df_profil_filtre_temp = df_tot_temp[df_tot_temp["z_temp"]>=-lim_depth].copy()
    # df_profil_filtre_Vs_rho = df_profil_filtre_Vs_rho[df_profil_filtre_Vs_rho["dt"].isin(temps_voulu)].copy()

    df_profil_filtre_temp = df_profil_filtre_temp[df_profil_filtre_temp["dt_temp"].isin(temps_voulu)].copy()
    # df_profil_filtre_temp["z_temp"] = df_profil_filtre_temp["z_temp"] - 4


        

    # Utilisation de pivot ----------------------------------------------------------------------------------------------------------------
    grid_temp = df_profil_filtre_temp.pivot(index='z_temp',columns='dt_temp',values='temp')
    grid_Vs = df_profil_filtre_Vs_rho.pivot(index='z',columns='dt',values='Vs')
    grid_rho = df_profil_filtre_Vs_rho.pivot(index='z',columns='dt',values='rho_vrai')
    # print("Grid shape:", grid.shape)

    # Plot --------------------------------------------------------------------------------------------------------------------------------
    plt.rcParams['font.size']=20
    fig, ax = plt.subplots(1,3, figsize=(20,5))

    im0 = ax[0].imshow(grid_temp, cmap="viridis",aspect='auto', origin='lower',extent=[0, 120,
                grid_temp.index.min(), grid_temp.index.max()])
    # ax[0].colorbar()
    ax[0].set_ylabel('Depth (m)')
    ax[0].set_xlabel('Day')
    ax[0].set_title('Temperature (°C)')

    im1 = ax[1].imshow(grid_Vs, cmap="viridis", aspect='auto', origin='lower')#,extent=[0, 120,            grid_Vs.index.min(), grid_Vs.index.max()])
    ax[1].set_yticklabels([])
    ax[1].set_xlabel('Day')
    ax[1].set_title(f'{Vp_Vs} velocity (m/s)')


    im2 = ax[2].imshow(grid_rho, cmap="viridis", aspect='auto', origin='lower',extent=[0, 120,100,250])#grid_rho.index.min(), grid_rho.index.max()])
    ax[2].set_yticklabels([])
    ax[2].set_xlabel('Day')
    ax[2].set_title('True resistivity '+ r"$(\Omega{}.m)$")

    fig.colorbar(im0, ax=ax[0])#, label="Temperature (°C)")
    cbar = fig.colorbar(im1, ax=ax[1])#, label="S wave velocity (m/s)")
    cbar.ax.ticklabel_format(style='sci', scilimits=(0,0))
    cbar2 = fig.colorbar(im2, ax=ax[2])#, label="True rho (Ohm.m)")
    cbar2.ax.ticklabel_format(style='sci', scilimits=(0,0))

    if len(barre_vertical) != 0 :
        for a in ax.flat:
            for i in barre_vertical:
                a.axvline(x=i, color="r", linestyle="--")
    plt.tight_layout()
    if path_save_fig :
        plt.savefig(path_save_fig)
    # plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def three_plot_observable_geophy_dharrma(path_result,debut_rp,jours_rp,pas_rp, facies, barre_vertical=[],path_save_fig =None,lim_depth = 2.1):
    '''
    Fonction qui permet de tracer les observable géophysique du modèle direct

    Entrée :
    path_result (str) : chemin absolue ammenant au dossier contenant les fichiers résultats à plotter.
    pluie (str) : Lettre permettant de savoir quelle scénario de pluie on veut plotter. Inutile quand diff = True.
    debut_rp (int) : Nombre (en jour) qui correspond au début de la représentation graphique
    jours_rp (int) : Nombre de jour de simulation que l'on veut plotter.
    pas_rp (int) : Temps (en jour) correspondant au pas de représentation
    list_mesure (list) : Liste permettant de mettre des barres verticales sur les graphs correspondant au mesure de profil.
    diff (bool) : Si diff = True la fonction plot la différence entre le scénario A et B. Sinon plot uniquement la pluie séléctionné avec "pluie"
    barre_vertical (list) : idem que liste mesure avec une représentation différente.
    path_save_fig(str) : Chemin permettant de sauvgarder la figure si l'on veut   
    '''
    

    path_temp = f'{path_result}/input_ginette/S_temperature_t.dat'
    path_PS_V = f'{path_result}/sismique/{facies}/output_SL_kk4_PS_v_Phase.dat'
    path_rho_app = f'{path_result}/elec/{facies}/rho_app_AB2.dat'

    temp_txt = np.loadtxt(path_temp)
    PS_V_txt = np.loadtxt(path_PS_V)
    rho_app_txt = np.loadtxt(path_rho_app)

    dt_temp = temp_txt[:,0]
    z_temp = temp_txt[:,1]
    temp = temp_txt[:,2]

    dt_PS_V = PS_V_txt[:,0]
    freq = PS_V_txt[:,1]
    PS_V = PS_V_txt[:,2]

    dt_rho = rho_app_txt[:,0]
    AB_2 = rho_app_txt[:,1]
    rho_app = rho_app_txt[:,2]

    temps_voulu = []

    df_tot_Vs = pd.DataFrame({'dt_PS_V': dt_PS_V,'freq': freq, 'PS_V': PS_V})
    df_tot_rho = pd.DataFrame({'dt_rho': dt_rho,'AB_2': AB_2,'rho_app' : rho_app})
    df_tot_temp = pd.DataFrame({'dt_temp': dt_temp,'z_temp':z_temp, 'temp': temp})

    for t_profil in range (debut_rp*86400,jours_rp*86400,pas_rp*86400):
            temps_voulu.append(t_profil)

    df_profil_filtre_Vs = df_tot_Vs[df_tot_Vs['freq']>=150].copy()
    df_profil_filtre_Vs = df_profil_filtre_Vs[df_profil_filtre_Vs['freq']<=250].copy()
    df_profil_filtre_rho = df_tot_rho.copy()
    df_profil_filtre_temp = df_tot_temp[df_tot_temp["z_temp"]>=-lim_depth].copy()
    # df_profil_filtre_Vs_rho = df_profil_filtre_Vs_rho[df_profil_filtre_Vs_rho["dt"].isin(temps_voulu)].copy()

    # Utilisation de pivot ----------------------------------------------------------------------------------------------------------------
    grid_temp = df_profil_filtre_temp.pivot(index='z_temp',columns='dt_temp',values='temp')
    grid_Vs = df_profil_filtre_Vs.pivot(index='freq',columns='dt_PS_V',values='PS_V')
    grid_rho = df_profil_filtre_rho.pivot(index='AB_2',columns='dt_rho',values='rho_app')
    # print("Grid shape:", grid.shape)

    # Plot --------------------------------------------------------------------------------------------------------------------------------
    plt.rcParams['font.size']=20
    fig, ax = plt.subplots(1,3, figsize=(20,5))

    im0 = ax[0].imshow(grid_temp, cmap="viridis",aspect='auto', origin='lower',extent=[0, 120,
            grid_temp.index.min(), grid_temp.index.max()])
    # ax[0].colorbar()
    ax[0].set_ylabel('Depth (m)')
    ax[0].set_xlabel('Day')
    ax[0].set_title('Temperature (°C)')

    im1 = ax[1].imshow(grid_Vs, cmap="viridis", aspect='auto', origin='lower',extent=[0, 120,grid_Vs.index.min(), grid_Vs.index.max()])
    # ax[1].set_yticklabels([])
    ax[1].set_xlabel('Day')
    ax[1].set_ylabel('Frequency (Hz)')
    ax[1].set_title('Phase velocity (m/s)')

    im2 = ax[2].imshow(grid_rho, cmap="viridis", aspect='auto', origin='lower',extent=[0, 120,
            grid_rho.index.min(), grid_rho.index.max()])
    # ax[2].set_yticklabels([])
    ax[2].set_xlabel('Day')
    ax[2].set_yscale('log')
    ax[2].set_ylabel('AB/2 (m)')
    ax[2].invert_yaxis()
    ax[2].set_title('Mesured resistivity ' + r"$(\Omega{}.m)$")

    fig.colorbar(im0, ax=ax[0])
    cbar = fig.colorbar(im1, ax=ax[1])
    # cbar.ax.ticklabel_format(style='sci', scilimits=(0,0))
    cbar2 = fig.colorbar(im2, ax=ax[2])
    # cbar2.ax.ticklabel_format(style='sci', scilimits=(0,0))

    if len(barre_vertical) != 0 :
        for a in ax.flat:
            for i in barre_vertical:
                a.axvline(x=i, color="r", linestyle="--")
    plt.tight_layout()
    if path_save_fig :
        plt.savefig(path_save_fig)
    # plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_profil_observable_dharrma(path_result,jour_profil,facies,representation = 1,path_fig = None,lim_depth=2.1):
    '''
    
    '''

    path_temp = f'{path_result}/input_ginette/S_temperature_t.dat'
    path_PS_V = f'{path_result}/sismique/{facies}/output_SL_kk4_PS_v_Phase.dat'
    path_rho_app = f'{path_result}/elec/{facies}/rho_app_AB2.dat'
    path_saturation = f'{path_result}/input_ginette/S_saturation_profil_t.dat'

    # path_temp = path_transient + f'/temp/{facies}/pluie_{pluie}/S_temperature_t.dat'
    # path_PS_V = path_transient + f'/sismique/pluie_{pluie}/{facies}/output_SL_kk4_PS_v_Phase.dat'
    # path_rho_app = path_transient + f'/elec/{facies}/pluie_{pluie}/rho_app_AB2.dat'
    # path_saturation = path_transient + f'/hydro/{facies}/pluie_{pluie}/S_saturation_profil_t.dat'


    color_map = copper(np.linspace(0, 1, len(jour_profil)))

    temp_txt = np.loadtxt(path_temp)
    PS_V_txt = np.loadtxt(path_PS_V)
    rho_app_txt = np.loadtxt(path_rho_app)
    sat_txt = np.loadtxt(path_saturation)

    # sat_static_txt = np.loadtxt(path_sat_static)
    # PS_V_static_txt = np.loadtxt(path_PS_V_static)
    # rho_app_static_txt = np.loadtxt(path_rho_app_static)

    # Hauteur WT --------------------------------------------------------------------------------------------------------------
    # hauteur_WT = np.loadtxt(path_WT).tolist()
    # hauteur_jours_WT = list()
    # for jour in jour_profil :
    #     h = hauteur_WT[jour*96]
    #     hauteur_jours_WT.append(-h+4) # Car le toit de la nappe est en hauteur pour ginette mes en profondeur pour le modèle
    

    dt_temp = temp_txt[:,0]
    z_temp = temp_txt[:,1]
    temp = temp_txt[:,2]

    dt_sat = sat_txt[:,0]
    z_sat = sat_txt[:,1]
    sat = sat_txt[:,2]

    dt_PS_V = PS_V_txt[:,0]
    freq = PS_V_txt[:,1]
    PS_V = PS_V_txt[:,2]

    dt_rho = rho_app_txt[:,0]
    AB_2 = rho_app_txt[:,1]
    rho_app = rho_app_txt[:,2]

    # hauteur_WT_sat = sat_static_txt[:,0]
    # z_sat_static = sat_static_txt[:,1]
    # sat_static = sat_static_txt[:,2]

    # hauteur_WT = PS_V_static_txt[:,0]
    # freq_static = PS_V_static_txt[:,1]
    # PS_V_static = PS_V_static_txt[:,2]

    # dt_rho_static = rho_app_static_txt[:,0]
    # AB_2_static = rho_app_static_txt[:,1]
    # rho_app_static = rho_app_static_txt[:,2]


    # temps_voulu = []

    df_tot_Vs = pd.DataFrame({'dt_PS_V': dt_PS_V,'freq': freq, 'PS_V': PS_V})
    df_tot_rho = pd.DataFrame({'dt_rho': dt_rho,'AB_2': AB_2,'rho_app' : rho_app})
    df_tot_temp = pd.DataFrame({'dt_temp': dt_temp,'z_temp':z_temp, 'temp': temp})
    df_tot_sat = pd.DataFrame({'dt_sat': dt_sat,'z_sat':z_sat, 'sat': sat})

    # df_tot_sat_static = pd.DataFrame({'hauteur_WT': hauteur_WT_sat,'z_sat_static': z_sat_static, 'sat_static': sat_static})
    # df_tot_Vs_static = pd.DataFrame({'hauteur_WT': hauteur_WT,'freq': freq_static, 'PS_V': PS_V_static})
    # df_tot_rho_static = pd.DataFrame({'dt_rho': dt_rho_static,'AB_2': AB_2_static,'rho_app' : rho_app_static})


    for i, t_profil in enumerate (jour_profil):
        color = color_map[i,:]
        t_jours_sec = t_profil*86400

        # Hauteur WT -------------------------------------------------------------------------------------------------
        # WT_profil = round(hauteur_jours_WT[i],4)
        # df_tot_Vs_static["hauteur_WT"] = df_tot_Vs_static["hauteur_WT"].round(4)
        # df_tot_sat_static["hauteur_WT"] = df_tot_sat_static["hauteur_WT"].round(4)

        # print(WT_profil)


        profil_temp = df_tot_temp[df_tot_temp["dt_temp"]==t_jours_sec]
        profil_sat = df_tot_sat[df_tot_sat["dt_sat"]==t_jours_sec]
        profil_rho = df_tot_rho[df_tot_rho["dt_rho"]==t_jours_sec]
        profil_Vs = df_tot_Vs[df_tot_Vs["dt_PS_V"]==t_jours_sec]

        
        # profil_sat_static = df_tot_sat_static[df_tot_sat_static["hauteur_WT"]==WT_profil]
        # profil_rho_static = df_tot_rho_static[df_tot_rho_static["dt_rho"]==t_jours_sec]
        # profil_Vs_static = df_tot_Vs_static[df_tot_Vs_static["hauteur_WT"]==WT_profil]

        
        if representation == 1 :
        # print(profil_temp)
            if i == 0 :
                plt.rcParams['font.size']=20
                fig,ax = plt.subplots(1,4,figsize=(15,5))
            
            ax[0].plot(profil_sat['sat'],profil_sat['z_sat'],color=color,label = f'Day {t_profil}')
            ax[0].set_ylabel('Depth (m)')
            # ax[0,0].set_xticklabels([])
            ax[0].legend(loc = 'lower right')
            ax[0].set_xlabel('Saturation (-)')

            ax[1].plot(profil_rho['rho_app'],profil_rho['AB_2'],color=color)
            ax[1].set_ylabel('AB/2 (m)')
            ax[1].set_yscale('log')
            ax[1].set_ylim(profil_rho['AB_2'].min(), profil_rho['AB_2'].max())
            ax[1].invert_yaxis()
            # ax[0,1].set_xticklabels([])
            ax[1].set_xlabel('Mesured resistivity ' + r"$(\Omega{}.m)$")


            ax[2].plot(profil_Vs['PS_V'],profil_Vs['freq'],color=color)
            ax[2].set_ylim(150,250)
            # ax[2].set_xlim(230,280)
            ax[2].set_ylabel('Fréquence (Hz)')
            # ax[0,2].set_xticklabels([])
            ax[2].set_xlabel ('Surface-wave velocity (m/s)')

            ax[3].plot(profil_temp['temp'],profil_temp['z_temp'],color=color)
            ax[3].set_ylabel('Depth (m)')
            # ax[0,2].set_xticklabels([])
            ax[3].set_xlabel ('Temperature (°C)')

            # ax[1,0].plot(profil_sat_static['sat_static'],profil_sat_static['z_sat_static'])
            # ax[1,0].set_ylabel('Depth (m)')
            # ax[1,0].set_xlabel('Saturation (-)')

            # ax[1,1].plot(profil_rho_static['rho_app'],profil_rho['AB_2'])
            # ax[1,1].set_ylabel('AB/2 (m)')
            # ax[1,1].set_yscale('log')
            # ax[1,1].set_ylim(profil_rho['AB_2'].min(), profil_rho['AB_2'].max())
            # ax[1,1].invert_yaxis()
            # ax[1,1].set_xlabel('Resistivity (Ohm.m)')

            # ax[1,2].plot(profil_Vs_static['PS_V'],profil_Vs_static['freq'])
            # ax[1,2].set_ylim(150,250)
            # ax[1,2].set_xlim(230,280)
            # ax[1,2].set_ylabel('Fréquence (Hz)')
            # ax[1,2].set_xlabel ('Velocity (m/s)')
        
        elif representation == 2 :
            if i == 0 :
                plt.rcParams['font.size']=15 
                fig,ax = plt.subplots(len(jour_profil),4,figsize=(15,10)) #15,10

                   
            if i == len(jour_profil)-1:
                ax[i,0].set_xlabel("Saturation (-)")
                ax[i,1].set_xlabel('Mesured resistivity ' + r"$(\Omega{}.m)$")
                ax[i,2].set_xlabel ('Surface-wave velocity (m/s)')
            else :
                ax[i,0].set_xticklabels([])
                ax[i,1].set_xticklabels([])
                ax[i,2].set_xticklabels([])
            
                
            ax[i,0].set_ylabel('Depth (m)')
            ax[i,1].set_ylabel('AB/2 (m)')
            ax[i,2].set_ylabel('Fréquence (Hz)')
            ax[i,3].set_ylabel('Depth (m)')


            # ax[i,0].set_title(f'Day {t_profil}')
            ax[i,0].plot(profil_sat['sat'],profil_sat['z_sat'],color = 'blue')
            # ax[i,0].plot(profil_sat_static['sat_static'],profil_sat_static['z_sat_static'], color = 'red',linestyle='--')
            ax[i,0].set_ylim(-lim_depth,0)
            # ax[i,0].set_xlabel("Saturation (-)")

            if i == 2 :
                ax[i,1].plot(profil_rho['rho_app'],profil_rho['AB_2'],color = 'blue',label = 'Transient')
                # ax[i,1].plot(profil_rho_static['rho_app'],profil_rho_static['AB_2'],color = 'red',linestyle='--', label = 'Satic')
            else :
                ax[i,1].plot(profil_rho['rho_app'],profil_rho['AB_2'],color = 'blue')
                # ax[i,1].plot(profil_rho_static['rho_app'],profil_rho_static['AB_2'],color = 'red',linestyle='--')
            ax[i,1].set_yscale('log')
            ax[i,1].set_ylim(profil_rho['AB_2'].min(), profil_rho['AB_2'].max())
            # if facies =='silt':
            #     ax[i,1].set_xlim(350,850)
            # elif facies == 'clay':
            #     ax[i,1].set_xlim(190,235)
            ax[i,1].invert_yaxis()
            # ax[i,1].set_xlabel('Resistivity (Ohm.m)')

            ax[i,2].plot(profil_Vs['PS_V'],profil_Vs['freq'],color = 'blue')
            # ax[i,2].plot(profil_Vs_static['PS_V'],profil_Vs_static['freq'],color = 'red',linestyle='--')
            ax[i,2].set_ylim(150,250)
            if facies =='silt':
                ax[i,2].set_xlim(230,270)
            elif facies == 'clay':
                ax[i,2].set_xlim(180,220)
            
            ax[i,3].plot(profil_temp['temp'],profil_temp['z_temp'],color = 'blue')
            ax[i,3].set_ylim(-lim_depth,0)
            ax2 = ax[i,3].twinx()
            ax2.set_yticklabels([])
            ax2.set_ylabel(f'Day {t_profil}')
            # ax[2,i].set_ylabel('Fréquence (Hz)')
            # ax[i,2].set_xlabel ('Velocity (m/s)')
                
    plt.tight_layout()
    if path_fig :
        plt.savefig(path_fig)
    # plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_profil_propriete_dharrma(path_result,jour_profil,facies,representation = 1,path_fig = None, Vp_Vs = 'Vs'):

    '''
    
    '''

    path_temp = f'{path_result}/input_ginette/S_temperature_t.dat'
    path_Vs = f'{path_result}/sismique/{facies}/output_SL_kk4_Vp_Vs.dat'
    path_rho_vrai = f'{path_result}/elec/{facies}/rho_vrai.dat'
    path_saturation = f'{path_result}/input_ginette/S_saturation_profil_t.dat'

        # path_sat_static = path_static + f'/{facies}/WT_scenario_{pluie}/hydro/Saturation.dat'
        # path_rho_vrai_static = path_static + f'/{facies}/WT_scenario_{pluie}/elec/rho_vrai_static_temperature.dat'
        # path_Vs_static = path_static + f'/{facies}/WT_scenario_{pluie}/sismique/output_SL_kk3_Vp_Vs.dat'

    temp_txt = np.loadtxt(path_temp)
    Vs_txt = np.loadtxt(path_Vs)
    rho_vrai_txt = np.loadtxt(path_rho_vrai)
    sat_txt = np.loadtxt(path_saturation)

    # sat_static_txt = np.loadtxt(path_sat_static)
    # Vs_static_txt = np.loadtxt(path_Vs_static)
    # rho_vrai_static_txt = np.loadtxt(path_rho_vrai_static)

    color_map = copper(np.linspace(0, 1, len(jour_profil)))

    # Hauteur WT --------------------------------------------------------------------------------------------------------------
    # hauteur_WT = np.loadtxt(path_WT).tolist()
    # hauteur_jours_WT = list()
    # for jour in jour_profil :
    #     h = hauteur_WT[jour*96]
    #     hauteur_jours_WT.append(-h+4) # Car le toit de la nappe est en hauteur pour ginette mes en profondeur pour le modèle
    

    dt_temp = temp_txt[:,0]
    z_temp = temp_txt[:,1]
    temp = temp_txt[:,2]

    dt_sat = sat_txt[:,0]
    z_sat = sat_txt[:,1]
    sat = sat_txt[:,2]

    dt_Vs = Vs_txt[:,0]
    z_sis = Vs_txt[:,1]
    if Vp_Vs == 'Vp':
        Vs = Vs_txt[:,2]
    else : 
        Vs = Vs_txt[:,3]

    dt_rho = rho_vrai_txt[:,0]
    z_elec = rho_vrai_txt[:,1]
    rho_vrai = rho_vrai_txt[:,2]

    # hauteur_WT_sat = sat_static_txt[:,0]
    # z_sat_static = sat_static_txt[:,1]
    # sat_static = sat_static_txt[:,2]

    # hauteur_WT = Vs_static_txt[:,0]
    # z_sis_static = Vs_static_txt[:,1]
    # Vs_static = Vs_static_txt[:,3]

    # dt_rho_static = rho_vrai_static_txt[:,0]
    # z_elec_static = rho_vrai_static_txt[:,1]
    # rho_vrai_static = rho_vrai_static_txt[:,2]

    

    # temps_voulu = []

    df_tot_Vs = pd.DataFrame({'dt_Vs': dt_Vs,'z_sis': z_sis, 'Vs': Vs})
    df_tot_rho = pd.DataFrame({'dt_rho': dt_rho,'z_elec': z_elec,'rho_vrai' : rho_vrai})
    df_tot_temp = pd.DataFrame({'dt_temp': dt_temp,'z_temp':z_temp, 'temp': temp})
    df_tot_sat = pd.DataFrame({'dt_sat': dt_sat,'z_sat':z_sat, 'sat': sat})

    # df_tot_sat_static = pd.DataFrame({'hauteur_WT': hauteur_WT_sat,'z_sat_static': z_sat_static, 'sat_static': sat_static})
    # df_tot_Vs_static = pd.DataFrame({'hauteur_WT': hauteur_WT,'z_sis_static': z_sis_static, 'Vs': Vs_static})
    # df_tot_rho_static = pd.DataFrame({'dt_rho': dt_rho_static,'z_elec_static': z_elec_static,'rho_vrai_static' : rho_vrai_static})


    for i, t_profil in enumerate (jour_profil):
        color = color_map[i,:]
        t_jours_sec = t_profil*86400

        # Hauteur WT -------------------------------------------------------------------------------------------------
        # WT_profil = round(hauteur_jours_WT[i],4)
        # df_tot_Vs_static["hauteur_WT"] = df_tot_Vs_static["hauteur_WT"].round(4)
        # df_tot_sat_static["hauteur_WT"] = df_tot_sat_static["hauteur_WT"].round(4)

        # print(WT_profil)


        profil_temp = df_tot_temp[df_tot_temp["dt_temp"]==t_jours_sec]
        profil_sat = df_tot_sat[df_tot_sat["dt_sat"]==t_jours_sec]
        profil_rho = df_tot_rho[df_tot_rho["dt_rho"]==t_jours_sec]
        profil_Vs = df_tot_Vs[df_tot_Vs["dt_Vs"]==t_jours_sec]

        
        # profil_sat_static = df_tot_sat_static[df_tot_sat_static["hauteur_WT"]==WT_profil]
        # profil_rho_static = df_tot_rho_static[df_tot_rho_static["dt_rho"]==t_jours_sec]
        # profil_Vs_static = df_tot_Vs_static[df_tot_Vs_static["hauteur_WT"]==WT_profil]

        
        if representation == 1 :
        # print(profil_temp)
            if i == 0 :
                plt.rcParams['font.size']=20
                fig,ax = plt.subplots(1,4,figsize=(15,5))
            
            ax[0].plot(profil_sat['sat'],profil_sat['z_sat'],color = color,label = f'Day {t_profil}')
            ax[0].set_ylabel('Depth (m)')
            # ax[0].legend(loc = 'lower right')
            # ax[0].set_xticklabels([])
            ax[0].set_xlabel('Saturation (-)')

            ax[1].plot(profil_rho['rho_vrai'],profil_rho['z_elec'],color = color)
            # ax[1].set_ylabel('AB/2 (m)')
            # ax[0,1].set_yscale('log')
            # ax[0,1].set_ylim(profil_rho['AB_2'].min(), profil_rho['AB_2'].max())
            # ax[1].invert_yaxis()
            # ax[1].set_xticklabels([])
            ax[1].set_xlabel('True Resistivity ' + r"$(\Omega{}.m)$")

            ax[2].plot(profil_Vs['Vs'],profil_Vs['z_sis'],color = color)
            ax[2].set_xlabel (f'{Vp_Vs} velocity (m/s)')
            # ax[2].set_ylim(150,250)
            # if facies == 'silt':
            #     ax[0,2].set_xlim(230,280)
            # # elif facies == 'clay':
            #     ax[0,2].set_xlim(230,275)
            # ax[2].set_ylabel('Fréquence (Hz)')
            # ax[2].set_xticklabels([])

            ax[3].plot(profil_temp['temp'],profil_temp['z_temp'],color=color)
            # ax[0,2].set_xticklabels([])
            ax[3].set_xlabel ('Temperature (°C)')

            # ax[1,0].plot(profil_sat_static['sat_static'],profil_sat_static['z_sat_static'])
            # ax[1,0].set_ylabel('Depth (m)')
            # ax[1,0].set_xlabel('Saturation (-)')

            # ax[1,1].plot(profil_rho_static['rho_app'],profil_rho['AB_2'])
            # ax[1,1].set_ylabel('AB/2 (m)')
            # ax[1,1].set_yscale('log')
            # ax[1,1].set_ylim(profil_rho['AB_2'].min(), profil_rho['AB_2'].max())
            # ax[1,1].invert_yaxis()
            # ax[1,1].set_xlabel('Resistivity ' + r"$(\Omega{}.m)$")

            # ax[1,2].plot(profil_Vs_static['PS_V'],profil_Vs_static['freq'])
            # ax[1,2].set_ylim(150,250)
            # ax[1,2].set_xlim(230,280)
            # ax[1,2].set_ylabel('Fréquence (Hz)')
            # ax[1,2].set_xlabel ('Velocity (m/s)')
        
        elif representation == 2 :
            if i == 0 :
                plt.rcParams['font.size']=20
                fig,ax = plt.subplots(len(jour_profil),4,figsize=(15,10)) #15,10
            ax[i,0].set_ylabel('Depth (m)')        

            ax[i,1].set_yticklabels([])
            ax[i,2].set_yticklabels([])
            ax[i,3].set_yticklabels([])
            
            if i == len(jour_profil)-1:
                ax[i,0].set_xlabel("Saturation (-)")
                ax[i,1].set_xlabel('True resistivity ' + r"$(\Omega{}.m)$")
                ax[i,2].set_xlabel ('S-wave velocity (m/s)')
                ax[i,3].set_xlabel('Temperature (°C)')

            else :
                ax[i,0].set_xticklabels([])
                ax[i,1].set_xticklabels([])
                ax[i,2].set_xticklabels([])
                ax[i,3].set_xticklabels([])

            
            
            ax[i,0].plot(profil_sat['sat'],profil_sat['z_sat'],color = 'blue')
            # ax[i,0].plot(profil_sat_static['sat_static'],profil_sat_static['z_sat_static'], color = 'red',linestyle='--')
            ax[i,0].set_ylim(-2.1,0)
            # ax[i,0].set_xlabel("Saturation (-)")

            if i == 1 :
                ax[i,1].plot(profil_rho['rho_vrai'],profil_rho['z_elec'],color = 'blue',label = 'Transient')
                # ax[i,1].plot(profil_rho_static['rho_vrai_static'],profil_rho_static['z_elec_static'],color = 'red',linestyle='--', label = 'Permanent')
                # ax[i,1].legend()
            else :
                ax[i,1].plot(profil_rho['rho_vrai'],profil_rho['z_elec'],color = 'blue')
                # ax[i,1].plot(profil_rho_static['rho_vrai_static'],profil_rho_static['z_elec_static'],color = 'red',linestyle='--')
            # ax[1,i].set_yscale('log')
            ax[i,1].set_ylim(-2.1,0)
            # if facies =='silt':
            #     ax[i,1].set_xlim(320,900)
            # elif facies =='clay':
            #     ax[i,1].set_xlim(185,245)
            # ax[1,i].invert_yaxis()
            # ax[i,1].set_xlabel('Resistivity (Ohm.m)')
                

            ax[i,2].plot(profil_Vs['Vs'],profil_Vs['z_sis'],color = 'blue')
            # ax[i,2].plot(profil_Vs_static['Vs'],profil_Vs_static['z_sis_static'],color = 'red',linestyle='--')
            ax[i,2].set_ylim(-2.1,0)
            if facies=='silt':
                ax[i,2].set_xlim(230,320)

            ax[i,3].plot(profil_temp['temp'],profil_temp['z_temp'],color = 'blue')
            ax2 = ax[i,3].twinx()
            ax2.set_yticklabels([])
            ax2.set_ylabel(f'Day {t_profil}')
            # ax[2,i].set_ylabel('Fréquence (Hz)')
            # ax[i,2].set_xlabel ('Velocity (m/s)')
                
    plt.tight_layout()
    if path_fig :
        plt.savefig(path_fig)
    # plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_profil_4_propriete(path_transient,pluie,path_static,path_WT,jour_profil,representation = 1,path_fig = None):

    '''
    
    '''
    path_temp = path_transient + f'/temp/pluie_{pluie}/S_temperature_t.dat'
    path_Vs = path_transient + f'/sismique/pluie_{pluie}/silt/output_SL_kk4_Vp_Vs.dat'
    path_rho_vrai = path_transient + f'/elec/pluie_{pluie}/rho_vrai.dat'
    path_saturation = path_transient + f'/hydro/pluie_{pluie}/S_saturation_profil_t.dat'
    path_vitesse = path_transient + f'/hydro/pluie_{pluie}//S_vitesse_profil.dat'

    path_sat_static = path_static + f'/WT_scenario_{pluie}/hydro/Saturation.dat'
    path_rho_vrai_static = path_static + f'/WT_scenario_{pluie}/elec/rho_vrai_static_temperature.dat'
    path_Vs_static = path_static + f'/WT_scenario_{pluie}/sismique/output_SL_kk3_Vp_Vs.dat'

    temp_txt = np.loadtxt(path_temp)
    Vs_txt = np.loadtxt(path_Vs)
    rho_vrai_txt = np.loadtxt(path_rho_vrai)
    sat_txt = np.loadtxt(path_saturation)

    sat_static_txt = np.loadtxt(path_sat_static)
    Vs_static_txt = np.loadtxt(path_Vs_static)
    rho_vrai_static_txt = np.loadtxt(path_rho_vrai_static)

    # Hauteur WT --------------------------------------------------------------------------------------------------------------
    hauteur_WT = np.loadtxt(path_WT).tolist()
    hauteur_jours_WT = list()
    for jour in jour_profil :
        h = hauteur_WT[jour*96]
        hauteur_jours_WT.append(-h+4) # Car le toit de la nappe est en hauteur pour ginette mes en profondeur pour le modèle
    

    dt_temp = temp_txt[:,0]
    z_temp = temp_txt[:,1]
    temp = temp_txt[:,2]

    dt_sat = sat_txt[:,0]
    z_sat = sat_txt[:,1]
    sat = sat_txt[:,2]

    dt_Vs = Vs_txt[:,0]
    z_sis = Vs_txt[:,1]
    Vs = Vs_txt[:,3]

    dt_rho = rho_vrai_txt[:,0]
    z_elec = rho_vrai_txt[:,1]
    rho_vrai = rho_vrai_txt[:,2]

    hauteur_WT_sat = sat_static_txt[:,0]
    z_sat_static = sat_static_txt[:,1]
    sat_static = sat_static_txt[:,2]

    hauteur_WT = Vs_static_txt[:,0]
    z_sis_static = Vs_static_txt[:,1]
    Vs_static = Vs_static_txt[:,3]

    dt_rho_static = rho_vrai_static_txt[:,0]
    z_elec_static = rho_vrai_static_txt[:,1]
    rho_vrai_static = rho_vrai_static_txt[:,2]

    

    # temps_voulu = []

    df_tot_Vs = pd.DataFrame({'dt_Vs': dt_Vs,'z_sis': z_sis, 'Vs': Vs})
    df_tot_rho = pd.DataFrame({'dt_rho': dt_rho,'z_elec': z_elec,'rho_vrai' : rho_vrai})
    df_tot_temp = pd.DataFrame({'dt_temp': dt_temp,'z_temp':z_temp, 'temp': temp})
    df_tot_sat = pd.DataFrame({'dt_sat': dt_sat,'z_sat':z_sat, 'sat': sat})

    df_tot_sat_static = pd.DataFrame({'hauteur_WT': hauteur_WT_sat,'z_sat_static': z_sat_static, 'sat_static': sat_static})
    df_tot_Vs_static = pd.DataFrame({'hauteur_WT': hauteur_WT,'z_sis_static': z_sis_static, 'Vs': Vs_static})
    df_tot_rho_static = pd.DataFrame({'dt_rho': dt_rho_static,'z_elec_static': z_elec_static,'rho_vrai_static' : rho_vrai_static})


    for i, t_profil in enumerate (jour_profil):
        t_jours_sec = t_profil*86400

        # Hauteur WT -------------------------------------------------------------------------------------------------
        WT_profil = round(hauteur_jours_WT[i],4)
        df_tot_Vs_static["hauteur_WT"] = df_tot_Vs_static["hauteur_WT"].round(4)
        df_tot_sat_static["hauteur_WT"] = df_tot_sat_static["hauteur_WT"].round(4)

        # print(WT_profil)


        profil_temp = df_tot_temp[df_tot_temp["dt_temp"]==t_jours_sec]
        profil_sat = df_tot_sat[df_tot_sat["dt_sat"]==t_jours_sec]
        profil_rho = df_tot_rho[df_tot_rho["dt_rho"]==t_jours_sec]
        profil_Vs = df_tot_Vs[df_tot_Vs["dt_Vs"]==t_jours_sec]

        
        profil_sat_static = df_tot_sat_static[df_tot_sat_static["hauteur_WT"]==WT_profil]
        profil_rho_static = df_tot_rho_static[df_tot_rho_static["dt_rho"]==t_jours_sec]
        profil_Vs_static = df_tot_Vs_static[df_tot_Vs_static["hauteur_WT"]==WT_profil]

        
        if representation == 1 :
        # print(profil_temp)
            if i == 0 :
                plt.rcParams['font.size']=20
                fig,ax = plt.subplots(2,3,figsize=(15,10))
            
            ax[0,0].plot(profil_sat['sat'],profil_sat['z_sat']-4)
            ax[0,0].set_ylabel('Depth (m)')
            ax[0,0].set_xticklabels([])

            ax[0,1].plot(profil_rho['rho_app'],profil_rho['AB_2'],label = f'Day {t_profil}')
            ax[0,1].set_ylabel('AB/2 (m)')
            # ax[0,1].set_yscale('log')
            ax[0,1].set_ylim(profil_rho['AB_2'].min(), profil_rho['AB_2'].max())
            ax[0,1].invert_yaxis()
            ax[0,1].set_xticklabels([])
            ax[0,1].legend(loc = 'lower right')

            ax[0,2].plot(profil_Vs['PS_V'],profil_Vs['freq'])
            # ax[0,2].set_ylim(150,250)
            # ax[0,2].set_xlim(230,280)
            ax[0,2].set_ylabel('Fréquence (Hz)')
            ax[0,2].set_xticklabels([])

            ax[1,0].plot(profil_sat_static['sat_static'],profil_sat_static['z_sat_static'])
            ax[1,0].set_ylabel('Depth (m)')
            ax[1,0].set_xlabel('Saturation (-)')

            ax[1,1].plot(profil_rho_static['rho_app'],profil_rho['AB_2'])
            ax[1,1].set_ylabel('AB/2 (m)')
            ax[1,1].set_yscale('log')
            ax[1,1].set_ylim(profil_rho['AB_2'].min(), profil_rho['AB_2'].max())
            ax[1,1].invert_yaxis()
            ax[1,1].set_xlabel('Resistivity (Ohm.m)')

            ax[1,2].plot(profil_Vs_static['PS_V'],profil_Vs_static['freq'])
            ax[1,2].set_ylim(150,250)
            ax[1,2].set_xlim(230,280)
            ax[1,2].set_ylabel('Fréquence (Hz)')
            ax[1,2].set_xlabel ('Velocity (m/s)')
        
        elif representation == 2 :
            if i == 0 :
                plt.rcParams['font.size']=20
                fig,ax = plt.subplots(len(jour_profil),3,figsize=(15,10)) #15,10
            ax[i,0].set_ylabel('Depth (m)')        

            ax[i,1].set_yticklabels([])
            ax[i,2].set_yticklabels([])
            
            if i == len(jour_profil)-1:
                ax[i,0].set_xlabel("Saturation (-)")
                ax[i,1].set_xlabel('True Resistivity (Ohm.m)')
                ax[i,2].set_xlabel ('Velocity S wave (m/s)')

            else :
                ax[i,0].set_xticklabels([])
                ax[i,1].set_xticklabels([])
                ax[i,2].set_xticklabels([])

            
            
            ax[i,0].plot(profil_sat['sat'],profil_sat['z_sat']-4,color = 'blue')
            ax[i,0].plot(profil_sat_static['sat_static'],profil_sat_static['z_sat_static'], color = 'red',linestyle='--')
            ax[i,0].set_ylim(-2.1,0)
            # ax[i,0].set_xlabel("Saturation (-)")

            if i == 1 :
                ax[i,1].plot(profil_rho['rho_vrai'],profil_rho['z_elec'],color = 'blue',label = 'Transient')
                ax[i,1].plot(profil_rho_static['rho_vrai_static'],profil_rho_static['z_elec_static'],color = 'red',linestyle='--', label = 'Permanent')
                # ax[i,1].legend()
            else :
                ax[i,1].plot(profil_rho['rho_vrai'],profil_rho['z_elec'],color = 'blue')
                ax[i,1].plot(profil_rho_static['rho_vrai_static'],profil_rho_static['z_elec_static'],color = 'red',linestyle='--')
            # ax[1,i].set_yscale('log')
            ax[i,1].set_ylim(-2.1,0)
            ax[i,1].set_xlim(320,900)
            # ax[1,i].invert_yaxis()
            # ax[i,1].set_xlabel('Resistivity (Ohm.m)')
                

            ax[i,2].plot(profil_Vs['Vs'],profil_Vs['z_sis'],color = 'blue')
            ax[i,2].plot(profil_Vs_static['Vs'],profil_Vs_static['z_sis_static'],color = 'red',linestyle='--')
            ax[i,2].set_ylim(-2.1,0)
            ax[i,2].set_xlim(230,320)
            ax2 = ax[i,2].twinx()
            ax2.set_yticklabels([])
            ax2.set_ylabel(f'Day {t_profil}')
            # ax[2,i].set_ylabel('Fréquence (Hz)')
            # ax[i,2].set_xlabel ('Velocity (m/s)')
                
    plt.tight_layout()
    if path_fig :
        plt.savefig(path_fig)
    plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_diff_relative_profil_observable(path_transient,pluie,path_static,path_WT,jour_profil,path_fig = None,facies='silt'):
    '''
    
    '''
    if facies =='silt':
        path_temp = path_transient + f'/temp/pluie_{pluie}/S_temperature_t.dat'
        path_PS_V = path_transient + f'/sismique/pluie_{pluie}/silt/output_SL_kk4_PS_v_Phase.dat'
        path_rho_app = path_transient + f'/elec/pluie_{pluie}/rho_app_AB2.dat'
        path_saturation = path_transient + f'/hydro/pluie_{pluie}/S_saturation_profil_t.dat'

        path_sat_static = path_static + f'/WT_scenario_{pluie}/hydro/Saturation.dat'
        path_rho_app_static = path_static + f'/WT_scenario_{pluie}/elec/rho_app_AB2_static_temperature.dat'
        path_PS_V_static = path_static + f'/WT_scenario_{pluie}/sismique/output_SL_kk3_PS_v_Phase.dat'
    
    else :
        path_temp = path_transient + f'/temp/{facies}/pluie_{pluie}/S_temperature_t.dat'
        path_PS_V = path_transient + f'/sismique/pluie_{pluie}/{facies}/output_SL_kk4_PS_v_Phase.dat'
        path_rho_app = path_transient + f'/elec/{facies}/pluie_{pluie}/rho_app_AB2.dat'
        path_saturation = path_transient + f'/hydro/{facies}/pluie_{pluie}/S_saturation_profil_t.dat'

        path_sat_static = path_static + f'/{facies}/WT_scenario_{pluie}/hydro/Saturation.dat'
        path_rho_app_static = path_static + f'/{facies}/WT_scenario_{pluie}/elec/rho_app_AB2_static_temperature.dat'
        path_PS_V_static = path_static + f'/{facies}/WT_scenario_{pluie}/sismique/output_SL_kk3_PS_v_Phase.dat'

    temp_txt = np.loadtxt(path_temp)
    PS_V_txt = np.loadtxt(path_PS_V)
    rho_app_txt = np.loadtxt(path_rho_app)
    sat_txt = np.loadtxt(path_saturation)

    sat_static_txt = np.loadtxt(path_sat_static)
    PS_V_static_txt = np.loadtxt(path_PS_V_static)
    rho_app_static_txt = np.loadtxt(path_rho_app_static)

    # Hauteur WT --------------------------------------------------------------------------------------------------------------
    hauteur_WT = np.loadtxt(path_WT).tolist()
    hauteur_jours_WT = list()
    for jour in jour_profil :
        h = hauteur_WT[jour*96]
        hauteur_jours_WT.append(-h+4) # Car le toit de la nappe est en hauteur pour ginette mes en profondeur pour le modèle
    

    dt_temp = temp_txt[:,0]
    z_temp = temp_txt[:,1]
    temp = temp_txt[:,2]

    dt_sat = sat_txt[:,0]
    z_sat = sat_txt[:,1]
    sat = sat_txt[:,2]

    dt_PS_V = PS_V_txt[:,0]
    freq = PS_V_txt[:,1]
    PS_V = PS_V_txt[:,2]

    dt_rho = rho_app_txt[:,0]
    AB_2 = rho_app_txt[:,1]
    rho_app = rho_app_txt[:,2]

    hauteur_WT_sat = sat_static_txt[:,0]
    z_sat_static = sat_static_txt[:,1]
    sat_static = sat_static_txt[:,2]

    hauteur_WT = PS_V_static_txt[:,0]
    freq_static = PS_V_static_txt[:,1]
    PS_V_static = PS_V_static_txt[:,2]

    dt_rho_static = rho_app_static_txt[:,0]
    AB_2_static = rho_app_static_txt[:,1]
    rho_app_static = rho_app_static_txt[:,2]

    # diff_sat = (sat-sat_static)/((sat+sat_static)/2)*100
    # (profil_sat['z_sat']-profil_sat_static['z_sat_static'])/((profil_sat['sat']-profil_sat_static['sat_static'])/2)*100

    # temps_voulu = []

    df_tot_Vs = pd.DataFrame({'dt_PS_V': dt_PS_V,'freq': freq, 'PS_V': PS_V})
    df_tot_rho = pd.DataFrame({'dt_rho': dt_rho,'AB_2': AB_2,'rho_app' : rho_app})
    df_tot_temp = pd.DataFrame({'dt_temp': dt_temp,'z_temp':z_temp, 'temp': temp})
    df_tot_sat = pd.DataFrame({'dt_sat': dt_sat,'z_sat':z_sat, 'sat': sat})

    df_tot_sat_static = pd.DataFrame({'hauteur_WT': hauteur_WT_sat,'z_sat_static': z_sat_static, 'sat_static': sat_static})
    df_tot_Vs_static = pd.DataFrame({'hauteur_WT': hauteur_WT,'freq': freq_static, 'PS_V': PS_V_static})
    df_tot_rho_static = pd.DataFrame({'dt_rho': dt_rho_static,'AB_2': AB_2_static,'rho_app' : rho_app_static})


    for i, t_profil in enumerate (jour_profil):
        t_jours_sec = t_profil*86400

        # Hauteur WT -------------------------------------------------------------------------------------------------
        WT_profil = round(hauteur_jours_WT[i],4)
        df_tot_Vs_static["hauteur_WT"] = df_tot_Vs_static["hauteur_WT"].round(4)
        df_tot_sat_static["hauteur_WT"] = df_tot_sat_static["hauteur_WT"].round(4)

        # print(WT_profil)


        profil_temp = df_tot_temp[df_tot_temp["dt_temp"]==t_jours_sec]
        profil_sat = df_tot_sat[df_tot_sat["dt_sat"]==t_jours_sec]
        profil_rho = df_tot_rho[df_tot_rho["dt_rho"]==t_jours_sec]
        profil_Vs = df_tot_Vs[df_tot_Vs["dt_PS_V"]==t_jours_sec]

        
        profil_sat_static = df_tot_sat_static[df_tot_sat_static["hauteur_WT"]==WT_profil]
        profil_rho_static = df_tot_rho_static[df_tot_rho_static["dt_rho"]==t_jours_sec]
        profil_Vs_static = df_tot_Vs_static[df_tot_Vs_static["hauteur_WT"]==WT_profil]

        

        
        if i == 0 :
                plt.rcParams['font.size']=20 
                fig,ax = plt.subplots(len(jour_profil),3,figsize=(15,10)) #15,10

                   
        if i == len(jour_profil)-1:
                ax[i,0].set_xlabel("Saturation (-)")
                ax[i,1].set_xlabel('Mesured resistivity ' + r"$(\Omega{}.m)$")
                ax[i,2].set_xlabel ('Surface-wave velocity (m/s)')
        # else :
        #         ax[i,0].set_xticklabels([])
        #         ax[i,1].set_xticklabels([])
        #         ax[i,2].set_xticklabels([])
            
                
        ax[i,0].set_ylabel('Depth (m)')
        ax[i,1].set_ylabel('AB/2 (m)')
        ax[i,2].set_ylabel('Fréquence (Hz)')

        # print (len(profil_sat['sat']),len(profil_sat_static['sat_static']),len(profil_sat['z_sat']))
        # diff_sat = (profil_sat['sat']-profil_sat_static['sat_static'])/((profil_sat['sat']+profil_sat_static['sat_static'])/2)*100
        # print(len(diff_sat))
        # print(diff_sat)
        # Interpolation de la saturation statique sur les profondeurs du profil dynamique
        sat_static_interp = np.interp(
            profil_sat['z_sat'], 
            profil_sat_static['z_sat_static'], 
            profil_sat_static['sat_static']
        )

        # Calcul de la différence relative (%)
        diff_sat = (profil_sat['sat'] - sat_static_interp) / ((profil_sat['sat'] + sat_static_interp) / 2) * 100
        ax[i,0].plot(diff_sat,profil_sat['z_sat']-4,color = 'blue')
        ax[i,0].set_ylim(-2.1,0)

        
        # Interpolation de la saturation statique sur les profondeurs du profil dynamique
        rho_static_interp = np.interp(
            profil_sat['z_sat'], 
            profil_sat_static['z_sat_static'], 
            profil_sat_static['sat_static']
        )

        # Calcul de la différence relative (%)
        diff_rho = (profil_rho['rho_app']-profil_rho_static['rho_app'])/((profil_rho['rho_app']+profil_rho_static['rho_app'])/2)*100#(profil_sat['sat'] - sat_static_interp) / ((profil_sat['sat'] + sat_static_interp) / 2) * 100
        ax[i,1].plot(diff_rho,profil_rho['AB_2'],color = 'blue')
        # ax[i,1].plot(profil_rho_static['rho_app'],profil_rho_static['AB_2'],color = 'red',linestyle='--')
        ax[i,1].set_yscale('log')
        ax[i,1].set_ylim(profil_rho['AB_2'].min(), profil_rho['AB_2'].max())
        ax[i,1].invert_yaxis()
            # ax[i,1].set_xlabel('Resistivity (Ohm.m)')

        # Interpolation du profil statique sur les profondeurs du profil dynamique
        Vs_static_interp = np.interp(
            profil_Vs['freq'],                # Profondeurs du profil dynamique
            profil_Vs_static['freq'],  # Profondeurs du profil statique
            profil_Vs_static['PS_V']          # Valeurs Vs statiques connues
        )

        # Calcul de la différence relative (%) correctement
        diff_Vs = (profil_Vs['PS_V'] - Vs_static_interp) / (
                    (profil_Vs['PS_V'] + Vs_static_interp) / 2
                ) * 100
        ax[i,2].plot(diff_Vs,profil_Vs['freq'],color = 'blue')
        # ax[i,2].plot(profil_Vs_static['PS_V'],profil_Vs_static['freq'],color = 'red',linestyle='--')
        ax[i,2].set_ylim(150,250)
        ax2 = ax[i,2].twinx()
        ax2.set_yticklabels([])
        ax2.set_ylabel(f'Day {t_profil}')
            # ax[2,i].set_ylabel('Fréquence (Hz)')
            # ax[i,2].set_xlabel ('Velocity (m/s)')
                
    plt.tight_layout()
    if path_fig :
        plt.savefig(path_fig)
    plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_temp(path_transient,pluie,jour_profil,path_fig=None,facies='silt'):
    if facies =='silt':
        path_temp = path_transient + f'/temp/pluie_{pluie}/S_temperature_t.dat'

    
    else :
        path_temp = path_transient + f'/temp/{facies}/pluie_{pluie}/S_temperature_t.dat'
        print(facies)
    temp_txt = np.loadtxt(path_temp)

    

    dt_temp = temp_txt[:,0]
    z_temp = temp_txt[:,1]
    temp = temp_txt[:,2]


    # temps_voulu = []
    df_tot_temp = pd.DataFrame({'dt_temp': dt_temp,'z_temp':z_temp, 'temp': temp})


    for i, t_profil in enumerate (jour_profil):
        t_jours_sec = t_profil*86400

        profil_temp = df_tot_temp[df_tot_temp["dt_temp"]==t_jours_sec]
        # print(profil_temp["temp"])
        print(profil_temp["temp"].iloc[0])

        if i == 0 :
                plt.rcParams['font.size']=20 
                fig,ax = plt.subplots(len(jour_profil),3,figsize=(15,10)) #15,10

                   
        # if i == len(jour_profil)-1:
        ax[i,0].set_xlabel("Temperature (°C)")
        # else :
        #         ax[i,0].set_xticklabels([])
        
        ax[i,0].set_ylabel('Depth (m)')

        ax[i,0].plot(profil_temp['temp'],profil_temp['z_temp']-4,color = 'blue',label = 'Transient')
        ax[i,0].plot([profil_temp['temp'].iloc[0],profil_temp['temp'].iloc[-1]],[profil_temp['z_temp'].iloc[0]-4,profil_temp['z_temp'].iloc[-1]-4],color = 'red',label = 'permanent',linestyle = '--')

    plt.tight_layout()
    if path_fig :
        plt.savefig(path_fig)
    plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_rho_Vp(path_transient,jours,pluie,facies):

    if facies =='silt':
        path_temp = path_transient + f'/temp/pluie_{pluie}/S_temperature_t.dat'
        path_Vp_Vs = path_transient + f'/sismique/pluie_{pluie}/silt/output_SL_kk4_Vp_Vs.dat'
        path_rho_vrai = path_transient + f'/elec/pluie_{pluie}/rho_vrai.dat'
        path_saturation = path_transient + f'/hydro/pluie_{pluie}/S_saturation_profil_t.dat'

    else :
        path_temp = path_transient + f'/temp/{facies}/pluie_{pluie}/S_temperature_t.dat'
        path_Vp_Vs = path_transient + f'/sismique/pluie_{pluie}/{facies}/output_SL_kk4_Vp_Vs.dat'
        path_rho_vrai = path_transient + f'/elec/{facies}/pluie_{pluie}/rho_vrai.dat'
        path_saturation = path_transient + f'/hydro/{facies}/pluie_{pluie}/S_saturation_profil_t.dat'
    
    temp_txt = np.loadtxt(path_temp)
    Vp_Vs_txt = np.loadtxt(path_Vp_Vs)
    rho_vrai_txt = np.loadtxt(path_rho_vrai)
    sat_txt = np.loadtxt(path_saturation)


    dt_temp = temp_txt[:,0]
    z_temp = temp_txt[:,1]
    temp = temp_txt[:,2]

    dt_sat = sat_txt[:,0]
    z_sat = sat_txt[:,1]
    sat = sat_txt[:,2]

    dt_rho = rho_vrai_txt[:,0]
    z_elec = rho_vrai_txt[:,1]
    rho_vrai = rho_vrai_txt[:,2]

    dt_sis = Vp_Vs_txt[:,0]
    z_sis = Vp_Vs_txt[:,1]
    Vp = Vp_Vs_txt[:,2]
    Vs = Vp_Vs_txt[:,3]

    df_tot_sis = pd.DataFrame({'dt_sis': dt_sis,'z_sis': z_sis, 'Vp': Vp, 'Vs': Vs})
    df_tot_rho = pd.DataFrame({'dt_rho': dt_rho,'z_elec': z_elec,'rho_vrai' : rho_vrai})

    profil_vp = df_tot_sis[df_tot_sis["dt_sis"]==jours*86400]
    profil_rho = df_tot_rho[df_tot_rho["dt_rho"]==jours*86400]

    fig,ax = plt.subplots(figsize=(15,10))
    ax.plot(profil_rho['rho_vrai'],profil_vp['Vp'],color = 'blue')
    ax.set_xlabel('Rho')
    ax.set_ylabel('Vp')
    ax.set_ylim(300,700)
    plt.show()

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_propriete_sismique(jour_profil,representation = 1,path_fig = None, facies ='silt'):

    '''
    
    '''
    path_Vs = 'sismique/silt/output_SL_kk4_Vp_Vs.dat'
    Vs_txt = np.loadtxt(path_Vs)
    dt_Vs = Vs_txt[:,0]
    z_sis = Vs_txt[:,1]
    Vs = Vs_txt[:,3]

    df_tot_Vs = pd.DataFrame({'dt_Vs': dt_Vs,'z_sis': z_sis, 'Vs': Vs})

    
    for i, t_profil in enumerate (jour_profil):
        t_jours_sec = t_profil*86400

        profil_Vs = df_tot_Vs[df_tot_Vs["dt_Vs"]==t_jours_sec]

        if i == 0:
            plt.rcParams['font.size']=20
            fig,ax = plt.subplots(len(jour_profil),3,figsize=(15,10)) #15,10
        ax[i,2].plot(profil_Vs['Vs'],profil_Vs['z_sis'],color = 'blue')

#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_only_sat(path_result,debut_rp,jours_rp,pas_rp,lim_depth = -2.1):
    '''
    
    '''

    sat_txt = np.loadtxt(f'{path_result}/input_ginette/S_saturation_profil_t.dat')
        
    dt = sat_txt[:,0]
    z = sat_txt[:,1]
    sat = sat_txt[:,2]

    temps_voulu = []

    df_tot_sat = pd.DataFrame({'dt_sat': dt,'z_sat':z, 'sat': sat})

    df_tot_sat = df_tot_sat[df_tot_sat["z_sat"]>=lim_depth].copy()

    sec_profil = [sec for sec in range(int(debut_rp),jours_rp,int(pas_rp))]
    color_map = copper(np.linspace(0, 1, len(sec_profil)))

    # df_tot_sat_static = pd.DataFrame({'hauteur_WT': hauteur_WT_sat,'z_sat_static': z_sat_static, 'sat_static': sat_static})
    # df_tot_Vs_static = pd.DataFrame({'hauteur_WT': hauteur_WT,'z_sis_static': z_sis_static, 'Vs': Vs_static})
    # df_tot_rho_static = pd.DataFrame({'dt_rho': dt_rho_static,'z_elec_static': z_elec_static,'rho_vrai_static' : rho_vrai_static})


    for i, t_profil in enumerate (sec_profil):
        color = color_map[i,:]

        profil_sat = df_tot_sat[df_tot_sat["dt_sat"]==t_profil]

        if i == 0 :
                plt.rcParams['font.size']=20
                fig,ax = plt.subplots(figsize=(5,15))
                # profil_sat = df_tot_sat[df_tot_sat["dt_sat"]==900]
                # ax.plot(profil_sat['sat'],profil_sat['z_sat'],color = color,label = f'Day {t_profil}')
            
        ax.plot(profil_sat['sat'],profil_sat['z_sat'],color = color,label = i)
        ax.set_ylabel('Depth (m)')
        ax.legend(loc = 'lower left')
        ax.set_xlabel('Saturation (-)')
    
#-------------------------------------------------------------------------------------------------------------------------------------------------------

def plot_temp_scenario(path_temp):
    temp = pd.read_csv(path_temp, header=None, sep='\s+', names=['T_top','T_bottom'])
    days = [i*900/86400 for i in range(0,len(temp['T_top']))]
    
    fig, ax = plt.subplots(1,figsize=(12, 4))
    fig.subplots_adjust(bottom = 0.258)
    fig.subplots_adjust(left = 0.08,top=0.96)
    ax.plot(days, temp['T_top'], label='Top Temperature', color='r', lw=1.5)
    ax.plot(days, temp['T_bottom'], label='Bottom Temperature', color='blue', ls='--', lw=1.5)
    ax.set_xlim(0,max(days))
    ax.set_ylim(min(min(temp['T_top']),min(temp['T_bottom']))-2,max(max(temp['T_top']),max(temp['T_bottom']))+2)
    ax.set_xlabel('Times (Days)', fontsize=16)
    ax.set_ylabel('Temperature\n(°C)', fontsize=16)
    fig.legend(fontsize = 16, loc='lower center', frameon=False, ncol=2)