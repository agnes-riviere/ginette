import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter
from matplotlib.cm import copper
import matplotlib.patches as patches
import matplotlib.colors as mcolors

##### Fonction support ###############

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

def creation_graph_err_sis(path_PS_V,t_err,pas=2):

    PS_V_txt = np.loadtxt(path_PS_V)

    dt_PS_V = PS_V_txt[:,0]
    freq = PS_V_txt[:,1]
    PS_V = PS_V_txt[:,2]

    df_tot_sur_wave = pd.DataFrame({'dt_PS_V': dt_PS_V,'freq': freq, 'PS_V': PS_V})

    df_sur_wave = df_tot_sur_wave[df_tot_sur_wave['dt_PS_V']==t_err]
    df_sur_wave = df_sur_wave[df_sur_wave['freq']<250]
    
    f_picked = df_sur_wave['freq'].to_numpy()
    v_picked = df_sur_wave['PS_V'].to_numpy()

    f_picked = f_picked[::pas]
    v_picked = v_picked[::pas]

    return(f_picked,lorentzian_error(v_picked,f_picked,0.25,96))

def creation_err_elec(eccartement,rho1,rho2,absolue = False):
    '''
    Valeur venant de la thèse de Gunther 2004 : 'Inversion Methods and Resolution Analysis for the 2D/3D Reconstruction of Resistivity Structures 
    from DC Measurements'
    '''
    err_list = list()
    err_coeff_al = pd.DataFrame({"AB_2": [1,2,3,4,5,7,9,101],"pourcentage":[4,3,2.6,2.5,2.4,2.3,2.2,2.1]})
    err_voltage = pd.DataFrame({"AB_2" : [5,12,101], 'pourcentage': [0.5,1,2]})

    if absolue :
        for i,AB in enumerate(eccartement):
            loc=0
            for k in err_coeff_al["AB_2"]:
                err = 0  
                if k >= AB :
                    err = err + (err_coeff_al['pourcentage'][loc]/100)*rho1[i]
                    err = err + (err_coeff_al['pourcentage'][loc]/100)*rho2[i]
                    break
                else :
                    loc += 1
            loc=0
            for j in err_voltage["AB_2"]:
                if j >= AB :
                    err = err + (err_voltage['pourcentage'][loc]/100)*rho1[i]
                    err = err + (err_voltage['pourcentage'][loc]/100)*rho2[i]
                    err_list.append(err)
                    break
                else :
                    loc += 1
    else : 
        for i,AB in enumerate(eccartement):
            loc=0
            for k in err_coeff_al["AB_2"]:
                err = 0  
                if k >= AB :
                    err = err + err_coeff_al['pourcentage'][loc]
                    break
                else :
                    loc += 1
            loc=0
            for j in err_voltage["AB_2"]:
                if j >= AB :
                    err = err + err_voltage['pourcentage'][loc]
                    err_list.append(err)
                    break
                else :
                    loc += 1
    return err_list

def err_elec(eccartement,rho):
    '''
    Valeur venant de la thèse de Gunther 2004 : 'Inversion Methods and Resolution Analysis for the 2D/3D Reconstruction of Resistivity Structures 
    from DC Measurements'
    '''
    err_list = list()
    err_coeff_al = pd.DataFrame({"AB_2": [1,2,3,4,5,7,9,101],"pourcentage":[4,3,2.6,2.5,2.4,2.3,2.2,2.1]})
    err_voltage = pd.DataFrame({"AB_2" : [5,12,101], 'pourcentage': [0.5,1,2]})

    for i,AB in enumerate(eccartement):
            loc=0
            for k in err_coeff_al["AB_2"]:
                err = 0  
                if k >= AB :
                    err = err + (err_coeff_al['pourcentage'][loc]/100)*rho[i]
                    break
                else :
                    loc += 1
            loc=0
            for j in err_voltage["AB_2"]:
                if j >= AB :
                    err = err + (err_voltage['pourcentage'][loc]/100)*rho[i]
                    err_list.append(err)
                    break
                else :
                    loc += 1
    return err_list


##### Fonction affichage ##############

def figure_scenario_infiltration(path_pluie_A,path_WT_9R,path_pluie_B,path_Evapo,path_no_Evapo,path_temp):

    plt.rcParams['font.size']=16

    dt_WT = [i/(24*4) for i in range(1,11520)]
    depth_WT_A = np.loadtxt(path_WT_9R)

    pluie_A = np.loadtxt(path_pluie_A)

    cumul_A_list = list()
    cumul_A = 0

    pluie_B = np.loadtxt(path_pluie_B)

    cumul_B_list = list()
    cumul_B = 0

    Evapo = np.loadtxt(path_Evapo)

    cumul_Evapo_list = list()
    cumul_Evapo = 0

    no_Evapo = np.loadtxt(path_no_Evapo)

    cumul_no_Evapo_list = list()
    cumul_no_Evapo = 0

    time = list()
    i=0

    for temps_sec in range(900,tot_day*86400,hydro_step):
            temps_jour = temps_sec/86400

            cumul_A = cumul_A + pluie_A[i]
            cumul_A_list.append(cumul_A)
            cumul_B = cumul_B + pluie_B[i]
            cumul_B_list.append(cumul_B)

            cumul_Evapo = cumul_Evapo + Evapo[i]
            cumul_Evapo_list.append(cumul_Evapo)
            cumul_no_Evapo = cumul_no_Evapo + no_Evapo[i]
            cumul_no_Evapo_list.append(cumul_no_Evapo)

            time.append(temps_jour)
            i=i+1
    
    temp = pd.read_csv(path_temp, header=None, sep='\s+', names=['T_top','T_bottom'])
    days = [i*900/86400 for i in range(0,len(temp['T_top']))]
    
    fig, ax = plt.subplots(6,figsize=(13, 13))
    fig.subplots_adjust(left = 0.094,top=0.976, bottom=0.124,right=0.889,hspace=0.29)
    ax[0].plot(days, temp['T_top'], label='Top Temperature', color='r', lw=1.5)
    ax[0].plot(days, temp['T_bottom'], label='Bottom Temperature', color='blue', ls='--', lw=1.5)
    ax[0].set_xlim(0,120)
    ax[0].set_ylim(-8,32)
    ax[0].set_xticklabels([])
    ax[0].set_ylabel('Temperature\n(°C)', fontsize=16)    


    ax[1].plot(time,[p*1e8 for p in pluie_A])
    ax[1].fill_between(time,[p*1e8 for p in pluie_A], 0, alpha=0.3, color="blue")
    ax[1].set_xlim(0,120)
    ax[1].set_ylim(0,4.3)
    ax[1].tick_params(axis='both', labelsize=16)
    ax[1].set_xticklabels([])
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[1].set_ylabel('Infiltration\n(1e-8 m/s)', fontsize = 16, color ='blue')
    ax[1].text(0.5, 0.80, '9R scenario', transform=ax[1].transAxes,fontsize=14, va='bottom', ha='center')

    ax2 = ax[1].twinx()
    ax2.plot(time,[c*10000 for c in cumul_A_list],color = 'black')
    ax2.set_ylim(0, 0.8)
    ax2.tick_params(axis='both', labelsize=16)
    ax2.set_ylabel('Rain\ncumulation\n(1e-4 m)', fontsize = 16)

    ax[2].plot(dt_WT,depth_WT_A,color ='blue')
    ax[2].set_xlim(0,120)
    ax[2].set_ylim(-2.01,-1.74)
    ax[2].set_xticklabels([])
    ax[2].set_ylabel('Water Table\ndepth (m)', fontsize = 16)

    ax[3].plot(time,[p*1e8 for p in pluie_B],color = '#800080')
    ax[3].fill_between(time, [p*1e8 for p in pluie_B], 0, alpha=0.3, color="#800080")
    ax[3].set_xlim(0,120)
    ax[3].set_ylim(0,4.3)# 3.2E-08)
    ax[3].tick_params(axis='both', labelsize=16)
    ax[3].set_xticklabels([])
    ax[3].set_ylabel('Infiltration\n(1e-8 m/s)', fontsize = 16,color = "purple")
    ax[3].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[3].text(0.5, 0.80, '2R scenario', transform=ax[3].transAxes,fontsize=14, va='bottom', ha='center')

    ax3 = ax[3].twinx()
    ax3.plot(time,[c*10000 for c in cumul_B_list],color = 'black')
    ax3.set_ylim(0, 0.8)
    ax3.tick_params(axis='both', labelsize=16)
    ax3.set_ylabel('Rain\ncumulation\n(1e-4 m)', fontsize = 16)

    ax[5].plot(time,[p*1e6 for p in Evapo],color="#B8D418C6")
    ax[5].fill_between(time,[p*1e6 for p in Evapo], 0, alpha=0.3, color="#B8D418C6")
    ax[5].set_xlim(0,120)
    ax[5].set_ylim(-0.2, 1.2)
    ax[5].tick_params(axis='both', labelsize=16)
    ax[5].set_xlabel('Time (day)',fontsize = 16)
    ax[5].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[5].set_ylabel('Effective\nrain\n(1e-6 m/s)', fontsize = 16, color ='#B8D418C6')
    ax[5].text(0.5, 0.80, 'Evapo scenario', transform=ax[5].transAxes,fontsize=14, va='bottom', ha='center')

    ax4 = ax[5].twinx()
    ax4.plot(time,[c*10000 for c in cumul_Evapo_list],color = 'black')
    ax4.tick_params(axis='both', labelsize=16)
    ax4.set_ylabel('Rain\ncumulation\n(1e-4 m)', fontsize = 16)

    ax[4].plot(time,[p*1e6 for p in no_Evapo],color = "#00802B")
    ax[4].fill_between(time, [p*1e6 for p in no_Evapo], 0, alpha=0.3, color="#00802B")
    ax[4].set_xlim(0,120)
    ax[4].set_ylim(0,1.2)
    ax[4].tick_params(axis='both', labelsize=16)
    ax[4].set_xticklabels([])
    ax[4].set_ylabel('Effective\nrain\n(1e-6 m/s)', fontsize = 16,color = "#00802B")
    ax[4].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax[4].text(0.5, 0.80, 'No evapo scenario', transform=ax[4].transAxes,fontsize=14, va='bottom', ha='center')

    ax5 = ax[4].twinx()
    ax5.plot(time,[c*10000 for c in cumul_no_Evapo_list],color = 'black')
    ax5.set_ylim(0, 0.8)
    ax5.tick_params(axis='both', labelsize=16)
    ax5.set_ylabel('Rain\ncumulation\n(1e-4 m)', fontsize = 16)

    for i in range(3):
        if i == 0 :
            ax[i].plot(1, -7.08, marker='v', color='g')
            ax[i].plot(30, -7.08, marker='v', color='g')        
            ax[i].plot(40, -7.08, marker='v', color='g')
            ax[i].plot(60, -7.08, marker='v', color='g')

        if i == 1 :
            ax[i].plot(1, 0.19, marker='v', color='g')
            ax[i].plot(30, 0.19, marker='v', color='g')        
            ax[i].plot(40, 0.19, marker='v', color='g')
            ax[i].plot(60, 0.19, marker='v', color='g')

        if i == 2:
            ax[i].plot(1, -2.0, marker='v', color='g')
            ax[i].plot(30, -2.0, marker='v', color='g')        
            ax[i].plot(40, -2.0, marker='v', color='g')
            ax[i].plot(60, -2.0, marker='v', color='g')
            
            ax[i].text(1, -2.05, 'D-1', ha='center', va='center',color = 'g')
            ax[i].text(30, -2.05, 'D-30', ha='center', va='center',color = 'g')
            ax[i].text(60, -2.05, 'D-60', ha='center', va='center',color = 'g')
            ax[i].text(40, -2.05, 'D-40', ha='center', va='center',color = 'g')

    rect = patches.Rectangle((25.5, 0.05),10,4.2,linewidth=1.5,edgecolor='purple',linestyle='--',facecolor='none')
    ax[1].add_patch(rect)
    rect = patches.Rectangle((25.5, -2.01),10,0.265,linewidth=1.5,edgecolor='purple',linestyle='--',facecolor='none')
    ax[2].add_patch(rect)
    rect = patches.Rectangle((25.5, 0.05),10,4.2,linewidth=1.5,edgecolor='purple',linestyle='--',facecolor='none')
    ax[3].add_patch(rect)

    rect_legend = patches.Patch(
    facecolor='none',
    edgecolor='purple', linestyle = '--',
    linewidth=1.5,
    label='Identical rain and\nwater level period'
    )
    all_handles = []
    all_labels  = []
    for axs in ax.flatten():     # axs = liste ou tableau des Axes
        h, l = axs.get_legend_handles_labels()
        all_handles += h
        all_labels  += l
        for child in axs.figure.axes:       # cherche tous les axes twins
            if child is not axs and child.bbox.bounds == axs.bbox.bounds:
                h2, l2 = child.get_legend_handles_labels()
                all_handles += h2
                all_labels  += l2
    triangle_J10 = Line2D([], [], marker='v', color='green',
                      linestyle='None', markersize=10)
    all_handles += [triangle_J10]
    all_labels  += ['Day of geophysical\nsurvey simulation']

    all_handles += [rect_legend]
    all_labels += ['Identical rain and\nwater level period']

    fig.legend(all_handles, all_labels,fontsize = 16,loc='lower center',framealpha=1,ncol=4,frameon=False)

def figure_prop_profil(path_result,path_static,path_WT,jour_profil):


 
    # path_temp = path_result+'9R/hydro/S_temperature_t.dat'
    path_Vs = path_result+'9R/seismic/output_SL_kk4_Vp_Vs.dat'
    path_rho_vrai = path_result+'9R/elec/rho_vrai.dat'
    path_saturation = path_result+'9R/hydro/S_saturation_profil_t.dat'
        

    path_sat_static2 =  path_static +'hydro/saturation_profil.dat'
    path_rho_vrai_static = path_static+'elec/rho_vrai_static_temperature.dat'
    path_Vs_static = path_static + 'seismic/output_SL_kk3_Vp_Vs.dat'

    # temp_txt = np.loadtxt(path_temp)
    Vs_txt = np.loadtxt(path_Vs)
    rho_vrai_txt = np.loadtxt(path_rho_vrai)
    sat_txt = np.loadtxt(path_saturation)


    sat_static_txt2 = np.loadtxt(path_sat_static2)
    Vs_static_txt = np.loadtxt(path_Vs_static)
    rho_vrai_static_txt = np.loadtxt(path_rho_vrai_static)

    # Hauteur WT --------------------------------------------------------------------------------------------------------------
    hauteur_WT = np.loadtxt(path_WT).tolist()
    hauteur_jours_WT = list()
    for jour in jour_profil :
        h = hauteur_WT[jour*96]
        hauteur_jours_WT.append(h)

    # dt_temp = temp_txt[:,0]
    # z_temp = temp_txt[:,1]
    # temp = temp_txt[:,2]

    dt_sat = sat_txt[:,0]
    z_sat = sat_txt[:,1]
    sat = sat_txt[:,2]

    dt_Vs = Vs_txt[:,0]
    z_sis = Vs_txt[:,1]
    Vs = Vs_txt[:,3]

    dt_rho = rho_vrai_txt[:,0]
    z_elec = rho_vrai_txt[:,1]
    rho_vrai = rho_vrai_txt[:,2]

    dt_sat_static = sat_static_txt2[:,0]
    z_sat_static2 = sat_static_txt2[:,1]
    sat_static2 = sat_static_txt2[:,2]

    hauteur_WT = Vs_static_txt[:,0]
    z_sis_static = Vs_static_txt[:,1]
    Vs_static = Vs_static_txt[:,3]

    dt_rho_static = rho_vrai_static_txt[:,0]
    z_elec_static = rho_vrai_static_txt[:,1]
    rho_vrai_static = rho_vrai_static_txt[:,2]

    

    # temps_voulu = []

    df_tot_Vs = pd.DataFrame({'dt_Vs': dt_Vs,'z_sis': z_sis, 'Vs': Vs})
    df_tot_rho = pd.DataFrame({'dt_rho': dt_rho,'z_elec': z_elec,'rho_vrai' : rho_vrai})
    # df_tot_temp = pd.DataFrame({'dt_temp': dt_temp,'z_temp':z_temp, 'temp': temp})
    df_tot_sat = pd.DataFrame({'dt_sat': dt_sat,'z_sat':z_sat, 'sat': sat})


    df_tot_sat_static2 = pd.DataFrame({'dt_sat_static': dt_sat_static,'z_sat_static': z_sat_static2, 'sat_static': sat_static2})
    df_tot_Vs_static = pd.DataFrame({'hauteur_WT': hauteur_WT,'z_sis_static': z_sis_static, 'Vs': Vs_static})
    df_tot_rho_static = pd.DataFrame({'dt_rho': dt_rho_static,'z_elec_static': z_elec_static,'rho_vrai_static' : rho_vrai_static})


    for i, t_profil in enumerate (jour_profil):
        t_jours_sec = t_profil*86400
        # print(t_jours_sec)

        # Hauteur WT -------------------------------------------------------------------------------------------------
        WT_profil = round(hauteur_jours_WT[i],4)
        df_tot_Vs_static["hauteur_WT"] = df_tot_Vs_static["hauteur_WT"].round(4)

        # profil_temp = df_tot_temp[df_tot_temp["dt_temp"]==t_jours_sec]
        profil_sat = df_tot_sat[df_tot_sat["dt_sat"]==t_jours_sec]
        profil_rho = df_tot_rho[df_tot_rho["dt_rho"]==t_jours_sec]
        profil_Vs = df_tot_Vs[df_tot_Vs["dt_Vs"]==t_jours_sec]
        profil_sat2 = df_tot_sat_static2[df_tot_sat_static2["dt_sat_static"]==t_jours_sec]

        profil_rho_static = df_tot_rho_static[df_tot_rho_static["dt_rho"]==t_jours_sec]
        profil_Vs_static = df_tot_Vs_static[df_tot_Vs_static["hauteur_WT"]==-WT_profil]
        
       
        if i == 0 :
            plt.rcParams['font.size']=16
            fig,ax = plt.subplots(len(jour_profil),3,figsize=(10,8)) #15,10
            fig.subplots_adjust(left=0.07,bottom=0.14,right=0.95,top=0.98)
        ax[i,0].set_ylabel('Depth (m)')        

        ax[i,1].set_yticklabels([])
        ax[i,2].set_yticklabels([])
            
        if i == len(jour_profil)-1:
            ax[i,0].set_xlabel("Saturation (-)")
            ax[i,2].set_xlabel('True resistivity ' + r"$(\Omega{}.m)$")
            ax[i,1].set_xlabel ('S-wave velocity (m/s)')

        else :
            ax[i,0].set_xticklabels([])
            ax[i,1].set_xticklabels([])
            ax[i,2].set_xticklabels([])

            
            
        ax[i,0].plot(profil_sat['sat'],profil_sat['z_sat'],color = 'blue')
        ax[i,0].plot(profil_sat2['sat_static'],profil_sat2['z_sat_static'], color = 'red',linestyle='--')
        ax[i,0].set_ylim(-2.1,0)

        if i == 1 :
            ax[i,2].plot(profil_rho['rho_vrai'],profil_rho['z_elec'],color = 'blue',label = 'Transient model')
            ax[i,2].plot(profil_rho_static['rho_vrai_static'],profil_rho_static['z_elec_static'],color = 'red',linestyle='--', label = 'Hydrostatic model')
        else :
            ax[i,2].plot(profil_rho['rho_vrai'],profil_rho['z_elec'],color = 'blue')
            ax[i,2].plot(profil_rho_static['rho_vrai_static'],profil_rho_static['z_elec_static'],color = 'red',linestyle='--')

        ax[i,2].set_ylim(-2.1,0)
        ax[i,1].plot(profil_Vs['Vs'],profil_Vs['z_sis'],color = 'blue')
        ax[i,1].plot(profil_Vs_static['Vs'],profil_Vs_static['z_sis_static'],color = 'red',linestyle='--')
        ax[i,1].set_ylim(-2.1,0)
        ax[i,1].set_xlim(230,320)

        ax2 = ax[i,2].twinx()
        ax2.set_yticklabels([])
        ax2.set_ylabel(f'Day {t_profil}')
    fig.legend(fontsize=16,loc='lower center',ncol=2,frameon=False,prop={'weight': 'normal'})

def figure_obs_profil(path_result,path_static,path_WT,jour_profil,errorbar = False):

    path_temp = path_result + '9R/hydro/S_temperature_t.dat'
    path_PS_V = path_result + '9R/seismic/output_SL_kk4_PS_v_Phase.dat'
    path_rho_app = path_result + '9R/elec/rho_app_AB2.dat'
    path_rho_app_static = path_static+'/elec/rho_app_AB2_static_temperature.dat'
    path_PS_V_static = path_static + 'seismic/output_SL_kk3_PS_v_Phase.dat'

    temp_txt = np.loadtxt(path_temp)
    PS_V_txt = np.loadtxt(path_PS_V)
    rho_app_txt = np.loadtxt(path_rho_app)

    PS_V_static_txt = np.loadtxt(path_PS_V_static)
    rho_app_static_txt = np.loadtxt(path_rho_app_static)

    # Hauteur WT --------------------------------------------------------------------------------------------------------------
    hauteur_WT = np.loadtxt(path_WT).tolist()
    hauteur_jours_WT = list()
    for jour in jour_profil :
        h = hauteur_WT[jour*96]
        hauteur_jours_WT.append(-h) # Valeur en hauteur négative pour correspondre à la convention de profondeur
    

    dt_temp = temp_txt[:,0]
    z_temp = temp_txt[:,1]
    temp = temp_txt[:,2]


    dt_PS_V = PS_V_txt[:,0]
    freq = PS_V_txt[:,1]
    PS_V = PS_V_txt[:,2]
    
    dt_rho = rho_app_txt[:,0]
    AB_2 = rho_app_txt[:,1]
    rho_app = rho_app_txt[:,2]


    hauteur_WT = PS_V_static_txt[:,0]
    freq_static = PS_V_static_txt[:,1]
    PS_V_static = PS_V_static_txt[:,2]

    dt_rho_static = rho_app_static_txt[:,0]
    AB_2_static = rho_app_static_txt[:,1]
    rho_app_static = rho_app_static_txt[:,2]


    # temps_voulu = []

    df_tot_Vs = pd.DataFrame({'dt_PS_V': dt_PS_V,'freq': freq, 'PS_V': PS_V})
    df_tot_rho = pd.DataFrame({'dt_rho': dt_rho,'AB_2': AB_2,'rho_app' : rho_app})
    df_tot_temp = pd.DataFrame({'dt_temp': dt_temp,'z_temp':z_temp, 'temp': temp})
    


    df_tot_Vs_static = pd.DataFrame({'hauteur_WT': hauteur_WT,'freq': freq_static, 'PS_V': PS_V_static})
    df_tot_rho_static = pd.DataFrame({'dt_rho': dt_rho_static,'AB_2': AB_2_static,'rho_app' : rho_app_static})


    for i, t_profil in enumerate (jour_profil):
        t_jours_sec = t_profil*86400

        # Hauteur WT -------------------------------------------------------------------------------------------------
        WT_profil = round(hauteur_jours_WT[i],4)
        df_tot_Vs_static["hauteur_WT"] = df_tot_Vs_static["hauteur_WT"].round(4)

        profil_temp = df_tot_temp[df_tot_temp["dt_temp"]==t_jours_sec]
        profil_rho = df_tot_rho[df_tot_rho["dt_rho"]==t_jours_sec]
        profil_Vs = df_tot_Vs[df_tot_Vs["dt_PS_V"]==t_jours_sec]

        profil_rho_static = df_tot_rho_static[df_tot_rho_static["dt_rho"]==t_jours_sec]
        profil_Vs_static = df_tot_Vs_static[df_tot_Vs_static["hauteur_WT"]==WT_profil]

        if i == 0 :
            plt.rcParams['font.size']=16 
            fig,ax = plt.subplots(len(jour_profil),3,figsize=(10,8))
            plt.subplots_adjust(wspace=0.65) 
            fig.subplots_adjust(left=0.072,bottom=0.166,right=0.957,top=0.983)

                   
        if i == len(jour_profil)-1:
            ax[i,2].set_xlabel('Measured resistivity \n' + r"$(\Omega{}.m)$")
            ax[i,1].set_xlabel ('Surface-wave \n velocity (m/s)')
            ax[i,0].set_xlabel("Temperature (°C)")
        else :
            ax[i,0].set_xticklabels([])
            ax[i,1].set_xticklabels([])
            ax[i,2].set_xticklabels([])
            
                
        ax[i,0].set_ylabel('Depth (m)')
        ax[i,2].set_ylabel('AB/2 (m)')
        ax[i,1].set_ylabel('Frequency\n(Hz)')

        if i == 2 :
            ax[i,2].plot(profil_rho['rho_app'],profil_rho['AB_2'],color = 'blue',label = 'Transient model')
            ax[i,2].plot(profil_rho_static['rho_app'],profil_rho_static['AB_2'],color = 'red',linestyle='--', label = 'Hydrostatic model')
        else :
            ax[i,2].plot(profil_rho['rho_app'],profil_rho['AB_2'],color = 'blue')
            ax[i,2].plot(profil_rho_static['rho_app'],profil_rho_static['AB_2'],color = 'red',linestyle='--')
        ax[i,2].set_yscale('log')
        ax[i,2].set_ylim(profil_rho['AB_2'].min(), profil_rho['AB_2'].max())
        ax[i,2].invert_yaxis()


        ax[i,1].plot(profil_Vs['PS_V'],profil_Vs['freq'],color = 'blue')
        ax[i,1].plot(profil_Vs_static['PS_V'],profil_Vs_static['freq'],color = 'red',linestyle='--')
        ax[i,1].set_ylim(0,250)

        ax[i,0].set_xlim(1,11.5)
        ax[i,1].set_xlim(225,300)
        ax[i,2].set_xlim(130,280)

        ax[i,0].plot(profil_temp['temp'],profil_temp['z_temp'],color = 'blue')
        ax[i,0].plot([profil_temp['temp'].iloc[0],profil_temp['temp'].iloc[-1]],[profil_temp['z_temp'].iloc[0],profil_temp['z_temp'].iloc[-1]],color = 'red',linestyle='--')
        ax[i,0].set_ylim(-2.1,0)
        ax2 = ax[i,2].twinx()
        ax2.set_yticklabels([])
        ax2.set_ylabel(f'Day {t_profil}')

        #Calcul et ajout de barre d'erreur
        if errorbar:

            x_temp_err = profil_temp['temp'].iloc[::15].to_numpy()
            y_temp_err = profil_temp['z_temp'].iloc[::15].to_numpy()

            ax[i,0].errorbar(x_temp_err,y_temp_err,xerr = 0.2,fmt='o',markersize=0,capsize =2,color = 'blue',alpha=0.5)


            x_tranient_err = profil_Vs['PS_V'].iloc[::25].to_numpy()
            y_tranient_err = profil_Vs['freq'].iloc[::25].to_numpy()
            val_tranient_err = lorentzian_error(x_tranient_err,y_tranient_err,0.25,96)

            ax[i,1].errorbar(x_tranient_err,y_tranient_err,xerr = val_tranient_err,fmt='o',markersize=0,capsize =2,color = 'blue',alpha=0.2)

            x_perm_err = profil_Vs_static['PS_V'].iloc[::25].to_numpy()
            y_perm_err = profil_Vs_static['freq'].iloc[::25].to_numpy()
            val_perm_err = lorentzian_error(x_perm_err,y_perm_err,0.25,96)

            ax[i,1].errorbar(x_perm_err,y_perm_err,xerr = val_perm_err,fmt='o',markersize=0,capsize =2,color = 'red',alpha=0.2)

            x_transient_err_elec = profil_rho['rho_app'].iloc[::5].to_numpy()
            y_transient_err_elec = profil_rho['AB_2'].iloc[::5].to_numpy()
            val_transient_err_elec = err_elec(y_transient_err_elec,x_transient_err_elec)

            ax[i,2].errorbar(x_transient_err_elec,y_transient_err_elec,xerr = val_transient_err_elec,fmt='o',markersize=0,capsize =2,color = 'blue',alpha=0.2)

            x_perm_err_elec = profil_rho_static['rho_app'].iloc[::5].to_numpy()
            y_perm_err_elec = profil_rho_static['AB_2'].iloc[::5].to_numpy()
            val_perm_err_elec = err_elec(y_perm_err_elec,x_perm_err_elec)

            ax[i,2].errorbar(x_perm_err_elec,y_perm_err_elec,xerr = val_perm_err_elec,fmt='o',markersize=0,capsize =2,color = 'red',alpha=0.2)

    fig.legend(fontsize=16,loc='lower center',ncol=2,frameon=False,prop={'weight': 'normal'})

def figure_comp_hydrological_state(path_WT_A,path_WT_B,path_sortie_ginette):

    dt_WT = [i/(24*4) for i in range(1,11520)]
    depth_WT_A = np.loadtxt(path_WT_A)
    depth_WT_B = np.loadtxt(path_WT_B)
    
    vit_eau_A = np.loadtxt(path_sortie_ginette+f'9R/hydro/S_vitesse_profil.dat')
    pression_A_txt = np.loadtxt(path_sortie_ginette+f'9R/hydro/S_pressure_profil_t.dat')
    sat_A_txt = np.loadtxt(path_sortie_ginette+f'9R/hydro/S_saturation_profil_t.dat')

    vit_eau_B = np.loadtxt(path_sortie_ginette+f'2R/hydro/S_vitesse_profil.dat')
    pression_B_txt = np.loadtxt(path_sortie_ginette+f'2R/hydro/S_pressure_profil_t.dat')
    sat_B_txt = np.loadtxt(path_sortie_ginette+f'2R/hydro/S_saturation_profil_t.dat')


    dt = vit_eau_A[:,0]
    z = vit_eau_A[:,1]
    vit_A = vit_eau_A[:,2]
    pression_A = pression_A_txt[:,2]
    sat_A = sat_A_txt[:,2]

    vit_B = vit_eau_B[:,2]
    pression_B = pression_B_txt[:,2]
    sat_B = sat_B_txt[:,2]

    temps_voulu = []

    df_tot = pd.DataFrame({'dt': dt,'z': z, 'vit_A': vit_A*1e8, 'pression_A': pression_A/10000, 'sat_A' : sat_A, 'vit_B': vit_B*1e8,
                            'pression_B': pression_B/10000, 'sat_B' : sat_B,'diff_vit': ((vit_A - vit_B))*1e8,
                            'diff_pression': ((pression_A - pression_B))/10000,
                            'diff_sat' : ((sat_A - sat_B))})
    
    for t_profil in range (1*86400,120*86400,1*86400):
            temps_voulu.append(t_profil)

    df_profil_filtre = df_tot[df_tot["z"]>=-2.1].copy()
    df_profil_filtre = df_profil_filtre[df_profil_filtre["dt"].isin(temps_voulu)].copy()

    grid_vit_A = df_profil_filtre.pivot(index='z',columns='dt',values='vit_A')
    grid_pression_A = df_profil_filtre.pivot(index='z',columns='dt',values='pression_A')
    grid_sat_A = df_profil_filtre.pivot(index='z',columns='dt',values='sat_A')

    grid_vit_B = df_profil_filtre.pivot(index='z',columns='dt',values='vit_B')
    grid_pression_B = df_profil_filtre.pivot(index='z',columns='dt',values='pression_B')
    grid_sat_B = df_profil_filtre.pivot(index='z',columns='dt',values='sat_B')
    
    grid_vit_diff = df_profil_filtre.pivot(index='z',columns='dt',values='diff_vit')
    grid_pression_diff = df_profil_filtre.pivot(index='z',columns='dt',values='diff_pression')
    grid_sat_diff = df_profil_filtre.pivot(index='z',columns='dt',values='diff_sat')

    # Plot --------------------------------------------------------------------------------------------------------------------------------
    plt.rcParams['font.size']=16
    plt.rcParams['axes.unicode_minus'] = False
    fig = plt.figure(figsize=(12,17))
    fig.subplots_adjust(left=0.078,bottom=0.12,right=0.9,top=0.944)
    gs = gridspec.GridSpec(5, 6, height_ratios=[1,0.0051,1,1,1],width_ratios=[1,1,0.05,0.6,1,0.05],hspace=0)

    #Haueur WT
    ax00 = fig.add_subplot(gs[0, 0])
    ax00.get_gridspec().update(hspace=0.5)
    ax00.plot(dt_WT,depth_WT_A,color ='blue')
    ax00.set_title('Scenario 9R',fontsize=16, pad = 20)
    ax00.set_ylabel('Depth (m)')
    ax00.set_xlabel('Time (Day)')
    ax00.set_ylim(-2.1,-1.8)
    ax00.set_xlim(0,120)
    ax00.text(-0.31, 1.15, 'a)', transform=ax00.transAxes,
        fontsize=14, va='bottom', ha='center')

    ax01 = fig.add_subplot(gs[0, 1])
    ax01.plot(dt_WT,depth_WT_B,color ='#800080')
    ax01.set_title('Scenario 2R',fontsize=16, pad = 20)
    ax01.set_xlabel('Time (Day)')
    ax01.set_yticklabels([])
    ax01.set_ylim(-2.1,-1.8)
    ax01.set_xlim(0,120)


    ### Resultat
    ax10 = fig.add_subplot(gs[2, 0])
    ax10.imshow(grid_vit_A, cmap="viridis",aspect='auto', origin='lower',extent=[0, 120,
                grid_vit_A.index.min(), grid_vit_A.index.max()])
    ax10.set_xticklabels([])
    ax10.set_ylabel('Depth(m)')
    ax10.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax10.text(-0.31, 1.15, 'b)', transform=ax10.transAxes,
        fontsize=14, va='bottom', ha='center')

    ax20 = fig.add_subplot(gs[3, 0])
    ax20.imshow(grid_pression_A, cmap="viridis", aspect='auto', origin='lower',extent=[0, 120,
            grid_pression_A.index.min(), grid_pression_A.index.max()])
    ax20.set_xticklabels([])
    ax20.set_ylabel('Depth (m)')
    ax20.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    ax30 = fig.add_subplot(gs[4, 0])
    ax30.imshow(grid_sat_A, cmap="viridis", aspect='auto', origin='lower',extent=[0, 120,
            grid_sat_A.index.min(), grid_sat_A.index.max()])
    ax30.set_xlabel('Time (Day)')
    ax30.set_ylabel('Depth (m)')
    ax30.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    

    ax11 = fig.add_subplot(gs[2, 1])
    im3 = ax11.imshow(grid_vit_B, cmap="viridis",aspect='auto', origin='lower',extent=[0, 120,
                grid_vit_A.index.min(), grid_vit_A.index.max()])
    ax11.set_xticklabels([])
    ax11.set_yticklabels([])


    ax21 = fig.add_subplot(gs[3, 1])
    im4 = ax21.imshow(grid_pression_B, cmap="viridis", aspect='auto', origin='lower',extent=[0, 120,
            grid_pression_A.index.min(), grid_pression_A.index.max()])
    ax21.set_xticklabels([])
    ax21.set_yticklabels([])

    ax31 = fig.add_subplot(gs[4, 1])
    im5 = ax31.imshow(grid_sat_B, cmap="viridis", aspect='auto', origin='lower',extent=[0, 120,
            grid_sat_A.index.min(), grid_sat_A.index.max()])
    ax31.set_yticklabels([])
    ax31.set_xlabel('Time (Day)',fontsize=16)

    ### Graph Diff
    ax13 = fig.add_subplot(gs[2, 4])
    ax13.set_title('Difference between\n scenarios',fontsize=16, pad=25)
    im6=ax13.imshow(grid_vit_diff, cmap="seismic",aspect='auto', origin='lower',vmin=-(max(abs(grid_vit_diff.values.min()),grid_vit_diff.values.max())),
                            vmax=(max(abs(grid_vit_diff.values.min()),grid_vit_diff.values.max())),extent=[0, 120,grid_vit_diff.index.min(), grid_vit_diff.index.max()])
    ax13.set_ylabel("Depth (m)")
    ax13.set_xticklabels([])
    ax13.text(-0.22, 1.15, 'c)', transform=ax13.transAxes,
        fontsize=14, va='bottom', ha='center')
    
    
    ax23 = fig.add_subplot(gs[3, 4])
    im7=ax23.imshow(grid_pression_diff, cmap="seismic",aspect='auto', origin='lower',vmin=-(max(abs(grid_pression_diff.values.min()),grid_pression_diff.values.max())),
                            vmax=(max(abs(grid_pression_diff.values.min()),grid_pression_diff.values.max())),extent=[0, 120,grid_pression_diff.index.min(), grid_pression_diff.index.max()])
    ax23.set_ylabel("Depth (m)")
    ax23.set_xticklabels([])
    
    ax33 = fig.add_subplot(gs[4, 4])
    im8=ax33.imshow(grid_sat_diff, cmap="seismic",aspect='auto', origin='lower',vmin=-(max(abs(grid_sat_diff.values.min()),grid_sat_diff.values.max())),
                            vmax=(max(abs(grid_sat_diff.values.min()),grid_sat_diff.values.max())),extent=[0, 120,grid_sat_diff.index.min(), grid_sat_diff.index.max()])
    ax33.set_ylabel("Depth (m)")
    ax33.set_xlabel('Time (Day)')

    cax12 = fig.add_subplot(gs[2, 2])
    cax22 = fig.add_subplot(gs[3, 2])
    cax32 = fig.add_subplot(gs[4, 2])
    cb = plt.colorbar(im3, cax=cax12,label='Water velocity\n 1e-8(m/s)')
    cb.ax.tick_params(labelsize=16)
    
    cb = plt.colorbar(im4, cax=cax22,label = "Pressure\n1e4(Pa)")
    cb.ax.tick_params(labelsize=16)
    cb.formatter = mticker.ScalarFormatter()
    cb.formatter.set_scientific(True)
    cb.formatter.set_powerlimits((0, 0))
    cb.update_ticks()

    cb = plt.colorbar(im5, cax=cax32,label ='Saturation (-)')
    cb.ax.tick_params(labelsize=16)
    

    cax14 = fig.add_subplot(gs[2, 5])
    cax24 = fig.add_subplot(gs[3, 5])
    cax34 = fig.add_subplot(gs[4, 5])
    cb = plt.colorbar(im6, cax=cax14,label='Diff velocity\n1e-8(m/s)')
    cb.ax.tick_params(labelsize=16)
    cb = plt.colorbar(im7, cax=cax24, label = "Diff pressure\n1e4(Pa)")
    cb.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    cb = plt.colorbar(im8, cax=cax34,label="Diff saturation\n(-)")
    cb.ax.tick_params(labelsize=16)

    axes = [ax00,ax01,ax10,ax11,ax13,ax20,ax21,ax23,ax30,ax31,ax33]
    
    for ax in axes :
        if ax == ax00 or ax == ax01:
            rect = patches.Rectangle((25.3, -2.1),10.7,0.296,linewidth=1.5,edgecolor='purple',linestyle='--',facecolor='none')
            ax.add_patch(rect)
        else:
            rect = patches.Rectangle((25.3, -2.09),10.7,2.06,linewidth=1.5,edgecolor='purple',linestyle='--',facecolor='none')
            ax.add_patch(rect)

    rect_legend = patches.Patch(
    facecolor='none',
    edgecolor='purple', linestyle = '--',
    linewidth=1.5,
    label='Identical rain and water level period'
    )

    fig.legend(handles=[rect_legend],loc='lower center',fontsize =16,frameon=False)

def figure_diff_obser(path_result):

    temp_A = np.loadtxt(path_result+'9R/hydro/S_temperature_t.dat')
    S_wave_velocity_A_txt = np.loadtxt(path_result+'9R/seismic/output_SL_kk4_PS_v_Phase.dat')
    elec_mesured_A_txt = np.loadtxt(path_result+'9R/elec/rho_app_AB2.dat')

    temp_B = np.loadtxt(path_result+'2R/hydro/S_temperature_t.dat')
    S_wave_velocity_B_txt = np.loadtxt(path_result+'2R/seismic/output_SL_kk4_PS_v_Phase.dat')
    elec_mesured_B_txt = np.loadtxt(path_result+'2R/elec/rho_app_AB2.dat')

    dt_temp = temp_A[:,0]
    z = temp_A[:,1]
    temperature_A = temp_A[:,2]
    dt_sis = S_wave_velocity_A_txt[:,0]
    freq = S_wave_velocity_A_txt[:,1]
    S_wave_velocity_A = S_wave_velocity_A_txt[:,2]
    dt_elec = elec_mesured_A_txt[:,0]
    AB_2 = elec_mesured_A_txt[:,1]
    elec_mesu_A = elec_mesured_A_txt[:,2]

    temperature_B = temp_B[:,2]
    S_wave_velocity_B = S_wave_velocity_B_txt[:,2]
    elec_mesu_B = elec_mesured_B_txt[:,2]

    df_temp = pd.DataFrame({'dt': dt_temp,'z': z,'diff_temp': temperature_A-temperature_B})
    df_elec = pd.DataFrame({'dt': dt_elec,'AB2':AB_2,'diff_elec':((elec_mesu_A-elec_mesu_B)),'elecA':elec_mesu_A,'elecB':elec_mesu_B})
    df_sis = pd.DataFrame({'dt': dt_sis,'freq': freq, 'diff_sis':((S_wave_velocity_A-S_wave_velocity_B))})

    
    df_temp_filtre = df_temp[df_temp["z"]>=-2.1].copy()
    df_sis_filtre1 = df_sis[df_sis["freq"]>=0].copy()
    df_sis_filtre = df_sis_filtre1[df_sis_filtre1["freq"]<=250].copy()


    grid_temp_diff = df_temp_filtre.pivot(index='z',columns='dt',values='diff_temp')
    grid_elec_diff = df_elec.pivot(index='AB2',columns='dt',values='diff_elec')
    grid_sis_diff = df_sis_filtre.pivot(index='freq',columns='dt',values='diff_sis')

    ### Calcul Err_sis
    f_picked,err_sis=creation_graph_err_sis(path_result+'9R/seismic/output_SL_kk4_PS_v_Phase.dat',86400*60)
    ### Calcul Err_ert
    df_elec_err = df_elec[df_elec['dt']==86400*119].copy() #Err calculé pour le 119ème jour (jour avec l'erreur la plus grande)
    err_ERT = creation_err_elec(df_elec_err['AB2'],df_elec_err['elecA'].values,df_elec_err['elecB'].values,absolue=True)

    # Plot --------------------------------------------------------------------------------------------------------------------------------
    plt.rcParams['font.size']=16
    fig = plt.figure(figsize=(18,5))
    
    fig.subplots_adjust(left=0.054,bottom=0.219,right=0.938,top=0.867) 
    gs_main = gridspec.GridSpec(1,3,figure=fig,wspace=0.6)

    axes = []

    for i in range(3):
        # Sous-grid dans chaque subplot principal
        gs_inner = gridspec.GridSpecFromSubplotSpec(1, 3,
            subplot_spec=gs_main[i],
            width_ratios=[0.3,1,0.05])
        ax0 = fig.add_subplot(gs_inner[0])
        ax1 = fig.add_subplot(gs_inner[1])
        ax2 = fig.add_subplot(gs_inner[2])

        axes.extend([ax0, ax1, ax2])
    
    axes[0].text(-1.1, 1.10, 'a)', transform=axes[0].transAxes,fontsize=20, va='bottom', ha='center')
    axes[0].plot([0.2,0.2],[-2.1,0])
    axes[0].set_ylim([-2.1,0])
    axes[0].set_xlim([0,0.3])
    axes[0].set_xticks(np.linspace(0, 0.2, 2))
    axes[0].set_ylabel('Depth (m)')
    axes[0].set_xlabel(r'$\varepsilon_{\mathrm{T}}$')

    v_extrem_t = max(abs(grid_temp_diff.values.min()),abs(grid_temp_diff.values.max()))
    
    im1=axes[1].imshow(grid_temp_diff, cmap="seismic",aspect='auto', origin='lower',extent=[0, 120,
                grid_temp_diff.index.min(), grid_temp_diff.index.max()],vmin = -v_extrem_t, vmax = v_extrem_t)
    axes[1].set_xlabel('Times (Day)')
    axes[1].set_yticklabels([])
    
    cbar = fig.colorbar(im1,cax=axes[2],label="Diff temperature (°C)")

    axes[3].text(-1.0, 1.10, 'b)', transform=axes[3].transAxes,fontsize=20, va='bottom', ha='center')
    axes[3].plot(err_sis,f_picked)
    axes[3].set_ylim([grid_sis_diff.index.min(),grid_sis_diff.index.max()])
    axes[3].set_xlim([0,20])
    axes[3].set_xticks(np.linspace(10, 10, 1))
    axes[3].set_ylabel('Frequency (Hz)')
    axes[3].set_xlabel(r'$\varepsilon_{\mathrm{V}}$')

    v_extrem_sis = max(abs(grid_sis_diff.values.min()),abs(grid_sis_diff.values.max()))
    im2 = axes[4].imshow(grid_sis_diff, cmap="seismic",aspect='auto', origin='lower',extent=[0, 120,grid_sis_diff.index.min(), grid_sis_diff.index.max()],vmin=-v_extrem_sis,vmax=v_extrem_sis)
    axes[4].set_xlabel('Times (Day)')
    axes[4].set_yticklabels([])
    

    cbar2 = fig.colorbar(im2,cax=axes[5],label='Diff surface-wave\nvelocity (m/s)')

    axes[6].text(-1.0, 1.10, 'c)', transform=axes[6].transAxes,fontsize=20, va='bottom', ha='center')
    axes[6].plot(err_ERT,df_elec_err['AB2'])
    axes[6].set_ylim([grid_elec_diff.index.min(),grid_elec_diff.index.max()])
    axes[6].set_xlim([0,25])
    axes[6].set_xticks(np.linspace(0, 20, 2))
    axes[6].set_ylabel('AB/2 (m)')
    axes[6].set_yscale('log')
    axes[6].invert_yaxis()
    axes[6].set_xlabel(r'$\varepsilon_{\mathrm{R}}$')

    v_elec = max(abs(grid_elec_diff.values.min()),abs(grid_elec_diff.values.max()))+2
    im3 = axes[7].imshow(grid_elec_diff, cmap="seismic",aspect='auto', origin='lower',extent=[0, 120,
                grid_elec_diff.index.min(), grid_elec_diff.index.max()],vmin=-v_elec,vmax=v_elec)
    axes[7].set_xlabel('Times (Day)')    
    axes[7].set_yscale('log')
    axes[7].set_ylim(df_elec['AB2'].min(), df_elec['AB2'].max())
    axes[7].invert_yaxis()
    axes[7].set_yticklabels([])

    cbar3 = fig.colorbar(im3,cax=axes[8],label='Diff measured resitivity\n'+r"$(\Omega{}.m)$")

    rect = patches.Rectangle((25.5, -2.09),10,2.07,linewidth=1.5,edgecolor='purple',linestyle='--',facecolor='none')
    axes[1].add_patch(rect)
    rect = patches.Rectangle((25.5, 16),10,232,linewidth=1.5,edgecolor='purple',linestyle='--',facecolor='none')
    axes[4].add_patch(rect)
    rect = patches.Rectangle((25.5, 1.55),10,99,linewidth=1.5,edgecolor='purple',linestyle='--',facecolor='none')
    axes[7].add_patch(rect)

    rect_legend = patches.Patch(
    facecolor='none',
    edgecolor='purple', linestyle = '--',
    linewidth=1.5,
    label='Identical rain and water level period'
    )

    fig.legend(handles=[rect_legend],loc="lower center",fontsize =16,frameon=False)

def figure_diff_evapo(path_result):

    temp_A = np.loadtxt(path_result+'Evapo/hydro/S_temperature_t.dat')
    S_wave_velocity_A_txt = np.loadtxt(path_result+'/Evapo/seismic/output_SL_kk4_PS_v_Phase.dat')
    elec_mesured_A_txt = np.loadtxt(path_result+'/Evapo/elec/rho_app_AB2.dat')

    temp_B = np.loadtxt(path_result+'/no_Evapo/hydro/S_temperature_t.dat')
    S_wave_velocity_B_txt = np.loadtxt(path_result+'/no_Evapo/seismic/output_SL_kk4_PS_v_Phase.dat')
    elec_mesured_B_txt = np.loadtxt(path_result+'/no_Evapo/elec/rho_app_AB2.dat')

    dt_temp = temp_A[:,0]
    z = temp_A[:,1]
    temperature_A = temp_A[:,2]
    dt_sis = S_wave_velocity_A_txt[:,0]
    freq = S_wave_velocity_A_txt[:,1]
    S_wave_velocity_A = S_wave_velocity_A_txt[:,2]
    dt_elec = elec_mesured_A_txt[:,0]
    AB_2 = elec_mesured_A_txt[:,1]
    elec_mesu_A = elec_mesured_A_txt[:,2]

    temperature_B = temp_B[:,2]
    S_wave_velocity_B = S_wave_velocity_B_txt[:,2]
    elec_mesu_B = elec_mesured_B_txt[:,2]


    df_temp = pd.DataFrame({'dt': dt_temp,'z': z,'diff_temp': temperature_A-temperature_B})
    df_elec = pd.DataFrame({'dt': dt_elec,'AB2':AB_2,'diff_elec':((elec_mesu_A-elec_mesu_B)),'elecA':elec_mesu_A,'elecB':elec_mesu_B})
    df_sis = pd.DataFrame({'dt': dt_sis,'freq': freq, 'diff_sis':((S_wave_velocity_A-S_wave_velocity_B))})

    
    df_temp_filtre = df_temp[df_temp["z"]>=-2.1].copy()
    df_sis_filtre1 = df_sis[df_sis["freq"]>=0].copy()
    df_sis_filtre = df_sis_filtre1[df_sis_filtre1["freq"]<=250].copy()


    grid_temp_diff = df_temp_filtre.pivot(index='z',columns='dt',values='diff_temp')
    grid_elec_diff = df_elec.pivot(index='AB2',columns='dt',values='diff_elec')
    grid_sis_diff = df_sis_filtre.pivot(index='freq',columns='dt',values='diff_sis')

    ### Calcul Err_sis
    f_picked,err_sis=creation_graph_err_sis(path_result+f'9R/seismic/output_SL_kk4_PS_v_Phase.dat',86400*60)

    ### Calcul Err_ert
    df_elec_err = df_elec[df_elec['dt']==86400*119].copy() #Err calculé pour le 119ème jour (jour avec l'eereur la plus grande)
    err_ERT = creation_err_elec(df_elec_err['AB2'],df_elec_err['elecA'].values,df_elec_err['elecB'].values,absolue=True)

    # Plot --------------------------------------------------------------------------------------------------------------------------------
    plt.rcParams['font.size']=16
    fig = plt.figure(figsize=(18,4))
    
    fig.subplots_adjust(left=0.06,bottom=0.19,right=0.917,top=0.85) 
    gs_main = gridspec.GridSpec(1,3,figure=fig,wspace=0.6)

    axes = []

    for i in range(3):
        # Sous-grid dans chaque subplot principal
        gs_inner = gridspec.GridSpecFromSubplotSpec(1, 3,
            subplot_spec=gs_main[i],
            width_ratios=[0.3,1,0.05])
        ax0 = fig.add_subplot(gs_inner[0])
        ax1 = fig.add_subplot(gs_inner[1])
        ax2 = fig.add_subplot(gs_inner[2])

        axes.extend([ax0, ax1, ax2])
    
    axes[0].text(-1.1, 1.10, 'a)', transform=axes[0].transAxes,fontsize=20, va='bottom', ha='center')
    axes[0].plot([0.2,0.2],[-2.1,0])
    axes[0].set_ylim([-2.1,0])
    axes[0].set_xlim([0,0.3])
    axes[0].set_xticks(np.linspace(0, 0.2, 2))
    axes[0].set_ylabel('Depth (m)')
    axes[0].set_xlabel(r'$\varepsilon_{\mathrm{T}}$')

    v_extrem_t = max(abs(grid_temp_diff.values.min()),abs(grid_temp_diff.values.max()))
    
    im1=axes[1].imshow(grid_temp_diff, cmap="seismic",aspect='auto', origin='lower',extent=[0, 120,
                grid_temp_diff.index.min(), grid_temp_diff.index.max()],vmin = -v_extrem_t, vmax = v_extrem_t)
    axes[1].set_xlabel('Times (Day)')
    axes[1].set_yticklabels([])
    
    cbar = fig.colorbar(im1,cax=axes[2],label="Diff temperature (°C)")

    axes[3].text(-1.0, 1.10, 'b)', transform=axes[3].transAxes,fontsize=20, va='bottom', ha='center')
    axes[3].plot(err_sis,f_picked)
    axes[3].set_ylim([grid_sis_diff.index.min(),grid_sis_diff.index.max()])
    axes[3].set_xlim([0,20])
    axes[3].set_xticks(np.linspace(10, 10, 1))
    axes[3].set_ylabel('Frequency (Hz)')
    axes[3].set_xlabel(r'$\varepsilon_{\mathrm{V}}$')

    v_extrem_sis = max(abs(grid_sis_diff.values.min()),abs(grid_sis_diff.values.max()))
    im2 = axes[4].imshow(grid_sis_diff, cmap="seismic",aspect='auto', origin='lower',extent=[0, 120,grid_sis_diff.index.min(), grid_sis_diff.index.max()],vmin=-v_extrem_sis,vmax=v_extrem_sis)
    axes[4].set_xlabel('Times (Day)')
    axes[4].set_yticklabels([])
    

    cbar2 = fig.colorbar(im2,cax=axes[5],label='Diff surface-wave\nvelocity (m/s)')

    axes[6].text(-1.0, 1.10, 'c)', transform=axes[6].transAxes,fontsize=20, va='bottom', ha='center') #-0.31, 1.15, 'b)'
    axes[6].plot(err_ERT,df_elec_err['AB2'])
    axes[6].set_ylim([grid_elec_diff.index.min(),grid_elec_diff.index.max()])
    axes[6].set_xlim([0,35])
    axes[6].set_xticks(np.linspace(0, 30, 2))
    axes[6].set_ylabel('AB/2 (m)')
    axes[6].set_yscale('log')
    axes[6].invert_yaxis()
    axes[6].set_xlabel(r'$\varepsilon_{\mathrm{R}}$')

    v_elec = max(abs(grid_elec_diff.values.min()),abs(grid_elec_diff.values.max()))+2
    im3 = axes[7].imshow(grid_elec_diff, cmap="seismic",aspect='auto', origin='lower',extent=[0, 120,
                grid_elec_diff.index.min(), grid_elec_diff.index.max()],vmin=-v_elec,vmax=v_elec)
    axes[7].set_xlabel('Times (Day)')    
    axes[7].set_yscale('log')
    axes[7].set_ylim(df_elec['AB2'].min(), df_elec['AB2'].max())
    axes[7].invert_yaxis()
    axes[7].set_yticklabels([])

    cbar3 = fig.colorbar(im3,cax=axes[8],label='Diff measured resitivity\n'+r"$(\Omega{}.m)$")

def figure_Annex_sismique(path_result,jours_voulu):

    temp = np.loadtxt(path_result+f'9R/hydro/S_temperature_t.dat')
    Surface_velocity_txt = np.loadtxt(path_result+'9R/seismic/output_SL_kk4_PS_v_Phase.dat')
    elec_mesured_txt = np.loadtxt(path_result+'9R/elec/rho_app_AB2.dat')

    saturation_txt = np.loadtxt(path_result+f'9R/hydro/S_saturation_profil_t.dat')
    elec_true_txt = np.loadtxt(path_result+f'9R/elec/rho_vrai.dat')
    Vp_Vs_txt = np.loadtxt(path_result+f'9R/seismic/output_SL_kk4_Vp_Vs.dat')
    firstArrP_firstArrS_txt = np.loadtxt(path_result+f'9R/seismic/output_SL_kk4_1AP_1AV_Phase.dat')

    dt_temp = temp[:,0]
    z = temp[:,1]
    temperature = temp[:,2]
    sat = saturation_txt[:,2]

    time_V = Vp_Vs_txt[:,0]
    zs = Vp_Vs_txt[:,1]
    Vp = Vp_Vs_txt[:,2]
    Vs = Vp_Vs_txt[:,3]
    bulk_density = Vp_Vs_txt[:,4]

    time_Arr = firstArrP_firstArrS_txt[:,0]
    offset = firstArrP_firstArrS_txt[:,1]
    firstArrP = firstArrP_firstArrS_txt[:,2]
    firstArrS = firstArrP_firstArrS_txt[:,3]
    dt_sis = Surface_velocity_txt[:,0]
    freq = Surface_velocity_txt[:,1]
    Surface_velocity = Surface_velocity_txt[:,2]

    dt_elec_true = elec_true_txt[:,0]
    z_elec = elec_true_txt[:,1]
    elec_true = elec_true_txt[:,2]
    dt_elec_mesured = elec_mesured_txt[:,0]
    AB_2 = elec_mesured_txt[:,1]
    elec_mesu = elec_mesured_txt[:,2]

    df_temp_obs = pd.DataFrame({'dt': dt_temp,'z': z,'temp': temperature,'saturation':sat})

    df_elec_true = pd.DataFrame({'dt': dt_elec_true,'z_elec':z_elec,'elec_true':elec_true})
    df_elec_obs = pd.DataFrame({'dt': dt_elec_mesured,'AB2':AB_2,'elec_mesu':elec_mesu})

    df_sis = pd.DataFrame({'dt': dt_sis,'freq': freq, 'surface_wave':Surface_velocity})
    df_vitesse = pd.DataFrame({'time':time_V,'zs':zs,'Vp':Vp,'Vs':Vs,'density':bulk_density})
    df_arrival = pd.DataFrame({'time':time_Arr,'offset':offset,'1AP':firstArrP,'1AS':firstArrS})


    df_temp_filtre = df_temp_obs[df_temp_obs["z"]>=-2.1].copy()

    df_elec_true_filtre = df_elec_true[df_elec_true["z_elec"]>=-2.1].copy()

    df_sis_filtre = df_sis[df_sis["freq"]<=250].copy()
    df_vitesse = df_vitesse[df_vitesse['zs']>=-2.1].copy()

    plt.rcParams['font.size']=16
    color_map = copper(np.linspace(0, 1, len(jours_voulu)))
    

    for i, day in enumerate (jours_voulu):
        color = color_map[i,:]
        day_seconds = day *86400

        df_temp_obs_day = df_temp_filtre[df_temp_filtre['dt'] == day_seconds].copy()
        df_elec_true_day= df_elec_true_filtre[df_elec_true_filtre['dt'] == day_seconds].copy()
        df_elec_obs_day = df_elec_obs[df_elec_obs['dt'] == day_seconds].copy()

        df_sis_day = df_sis_filtre[df_sis_filtre['dt'] == day_seconds].copy()
        df_vitesse_day = df_vitesse[df_vitesse['time'] == day_seconds].copy()
        df_arrival_day = df_arrival[df_arrival['time'] == day_seconds].copy()

        if i == 0:
            fig, ax = plt.subplots(2,5,figsize = (18,10))
            fig.subplots_adjust(left=0.057,bottom=0.117,right=0.978,top=0.947,wspace=0.426,hspace=0.3)
        
        ax[0,0].plot(df_temp_obs_day['saturation'].values,df_temp_obs_day['z'].values,color=color,label = f'Day {day}')
        ax[0,0].set_ylabel('Depth (m)')
        ax[0,0].set_xlabel('Saturation (-)')
        ax[0,0].set_ylim(-2.1,0.0)
        ax[0,0].text(-0.35, 1.07, 'a)', transform=ax[0,0].transAxes,fontsize=16, va='bottom', ha='center')

        ax[0,1].plot(df_elec_true_day['elec_true'].values,df_elec_true_day['z_elec'].values,color=color)
        ax[0,1].set_xlabel('True resistivity ' + r"$(\Omega{}.m)$")
        ax[0,1].set_ylabel('Depth (m)')
        ax[0,1].set_ylim(-2.1,0.0)
        
        ax[0,2].plot(df_vitesse_day['Vp'].values,df_vitesse_day['zs'].values, color=color)
        ax[0,2].set_xlabel('P-waves velocity (m/s)')
        ax[0,2].set_ylabel('Depth (m)')
        ax[0,2].set_ylim(-2.1,0.0)

        ax[0,3].plot(df_vitesse_day['Vs'].values,df_vitesse_day['zs'].values, color=color)
        ax[0,3].set_xlabel('S-wave velocity (m/s)')
        ax[0,3].set_ylabel('Depth (m)')
        ax[0,3].set_ylim(-2.1,0.0)

        ax[0,4].plot(df_vitesse_day['density'].values,df_vitesse_day['zs'].values, color=color)
        ax[0,4].set_xlabel('Bulk density (Pa)')
        ax[0,4].set_ylabel('Depth (m)')
        ax[0,4].set_ylim(-2.1,0.0)

        ax[1,0].plot(df_temp_obs_day['temp'].values,df_temp_obs_day['z'].values,color=color)
        ax[1,0].set_ylabel('Depth (m)')
        ax[1,0].set_xlabel('Temperature (°C)')
        ax[1,0].set_ylim(-2.1,0.0)
        ax[1,0].text(-0.35, 1.07, 'b)', transform=ax[1,0].transAxes,
        fontsize=16, va='bottom', ha='center')

        ax[1,1].plot(df_elec_obs_day['elec_mesu'].values,df_elec_obs_day['AB2'].values,color=color)
        ax[1,1].set_ylim(df_elec_obs_day['AB2'].values.min(), df_elec_obs_day['AB2'].values.max())
        ax[1,1].invert_yaxis()
        ax[1,1].set_ylabel("AB/2")
        ax[1,1].set_xlabel("Measured resistivity " + r"$(\Omega{}.m)$")
        ax[1,1].set_yscale('log')

        ax[1,2].plot(df_arrival_day['offset'].values,df_arrival_day['1AP'].values, color=color)
        ax[1,2].set_ylabel('First P-wave arrival time (s)')
        ax[1,2].set_xlabel('Offset (m)')
        ax[1,2].set_ylim(0,df_arrival_day['1AS'].values.max())
        ax[1,2].set_xlim(0,df_arrival_day['offset'].values.max())

        ax[1,3].plot(df_arrival_day['offset'].values,df_arrival_day['1AS'].values, color=color)
        ax[1,3].set_ylabel('First S-wave arrival time (s)')
        ax[1,3].set_xlabel('Offset (m)')
        ax[1,3].set_ylim(0,df_arrival_day['1AS'].values.max())
        ax[1,3].set_xlim(0,df_arrival_day['offset'].values.max())

        ax[1,4].plot(df_sis_day['surface_wave'].values,df_sis_day['freq'].values, color=color)
        ax[1,4].set_ylabel("Frequency (Hz)")
        ax[1,4].set_xlabel("Surface wave velocity (m/s)")
        ax[1,4].set_ylim(0,250)
        ax[1,4].set_xlim(225,300)

    fig.legend(loc='lower center',ncol=4,fontsize =16,frameon=False)

###### Parametre ######################

tot_day = 120
hydro_step = 900

jour_profil = [1,30,40,60]


path_temp = 'data/E_temp_t.dat'
path_infiltration_9R ='data/transient/9R/hydro/E_debit_haut_t.dat'
path_infiltration_2R = 'data/transient/2R/hydro/E_debit_haut_t.dat'
path_infiltration_Evapo = 'data/transient/Evapo/hydro/E_debit_haut_t.dat'
path_infiltration_no_Evapo = 'data/transient/no_Evapo/hydro/E_debit_haut_t.dat'

path_data = 'data/transient/'
path_static ='data/hydrostatic/'

path_WT_level_9R = 'data/transient/9R/hydro/Hauteur_WT_9R.dat'
path_WT_level_2R = 'data/transient/2R/hydro/Hauteur_WT_2R.dat'


###### Lancement fonction ##############

# Figure 2
figure_scenario_infiltration(path_infiltration_9R,path_WT_level_9R,path_infiltration_2R,path_infiltration_Evapo,
                             path_infiltration_no_Evapo,path_temp)

# Figure 4
figure_prop_profil(path_data,path_static,path_WT_level_9R,jour_profil)

# Figure 5
figure_obs_profil(path_data,path_static,path_WT_level_9R,jour_profil,errorbar=True)

# Figure 6
figure_comp_hydrological_state(path_WT_level_9R,path_WT_level_2R,path_data)

# Figure 7
figure_diff_obser(path_data)

# Figure 8
figure_diff_evapo(path_data)

# Annexe seismic
figure_Annex_sismique(path_data,jour_profil)


plt.show()