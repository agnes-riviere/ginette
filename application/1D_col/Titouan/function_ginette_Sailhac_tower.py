import os
import pandas as pd
import numpy as np
# from IPython.display import display

os.chdir('/home/th/Documents/code/ginette/application/1D_col')
# os.chdir('/home/el/Codes/Titouan/ginette/application/1D_col')
         
def preprocess(nb_day,dt,nb_zone):
    ########### Setup of the model

    #time step in s
    dt=dt
    #duration of the simulation in days
    nb_day=nb_day

    # state
    ## 0 steady state
    # 1 transient state (dynamic state)
    state=1

    # altitude column top in meter
    z_top=0
    #altitude column bottom in meter
    z_bottom=-1.8
    az=abs(z_top-z_bottom)
    #discretisation : size cell in meter
    dz=0.01
    
    #Observation in meter
    Obs1=z_top-0.025
    Obs2=z_top-0.05
    Obs3=z_top-0.1
    Obs4=z_top-0.2
    Obs5=z_top-0.3
    Obs6=z_top-0.5
    Obs7=z_top-0.7
    Obs8=z_top-1
    Obs9=z_top-1.5
    Obs10=z_bottom
    
    # thk2=prof_zone2

    #apply unsaturated flow and thermal 
    #unsat =1 apply
    #unsat=0 cancel unsaturated zone
    unsat=1

    # number of facies in the column. If nb_zone=1 homognous porous media
    nb_zone=nb_zone

    # number of cell
    ## write the parameters
    nb_cell=abs(z_top-z_bottom)/dz
    cell1=-Obs1/dz
    cell2=-Obs2/dz
    cell3=-Obs3/dz
    cell4=-Obs4/dz
    cell5=-Obs5/dz
    cell6=-Obs6/dz
    cell7=-Obs7/dz
    cell8=-Obs8/dz
    cell9=-Obs9/dz
    cell10=-Obs10/dz
    #-----------------------------------------------------------------
    ## write the setup of the moddeled domain    
    f_param_bck=open("E_parametre_backup.dat", "r")
    f_param_new = open("E_parametre.dat", 'w')
    setup_model=f_param_bck.read()
    setup_model=setup_model.replace('[dt]','%06.0fD+00' % dt)
    setup_model=setup_model.replace('[state]','%1i' % state)
    setup_model=setup_model.replace('[nb_day]','%06.0f' % nb_day)
    setup_model=setup_model.replace('[z_top]', '%7.3e' % z_top)
    setup_model=setup_model.replace('[z_bottom]', '%7.3e' % z_bottom)
    setup_model=setup_model.replace('[az]','%7.3e' % az)
    setup_model=setup_model.replace('[dz]','%6.2e' % dz)
    setup_model=setup_model.replace('[nb_cell]','%05.0f' % nb_cell)
    setup_model=setup_model.replace('[unsat]','%1i' % unsat)
    setup_model=setup_model.replace('[cell1]','%05d' % cell1)
    setup_model=setup_model.replace('[cell2]','%05d' % cell2)
    setup_model=setup_model.replace('[cell3]','%05d' % cell3)
    setup_model=setup_model.replace('[cell4]','%05d' % cell4)
    setup_model=setup_model.replace('[cell5]','%05d' % cell5)
    setup_model=setup_model.replace('[cell6]','%05d' % cell6)
    setup_model=setup_model.replace('[cell7]','%05d' % cell7)
    setup_model=setup_model.replace('[cell8]','%05d' % cell8)
    setup_model=setup_model.replace('[cell9]','%05d' % cell9)
    setup_model=setup_model.replace('[cell10]','%05d' % cell10)
    
    f_param_new.write(setup_model)
    f_param_bck.close()
    f_param_new.close()
    
        ########### Zone of parameters
    f_coor=open("E_coordonnee.dat", "w")
    f_zone=open("E_zone.dat", 'w')
    coord=pd.DataFrame()
    # calculate the coordinates to know the number of cell in the domain    
    zvalues = np.arange(-dz/2, z_bottom,-dz );
    xvalues = np.array([0.5]);
    zz, xx = np.meshgrid(zvalues, xvalues)
    NT = np.product(zz.shape)
    data = {
        "x": np.reshape(xx,NT),
        "z": np.reshape(zz,NT)}
    coord = pd.DataFrame(data=data)
    coord['id']=coord.index.values.astype(int)
    coord['id']=coord['id']+1
    cols = coord.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    coord = coord[cols] 
    coord.to_csv(f_coor, index = False, sep=' ', header=False)
    #zone parameter by cell ((homogenous domain = 1 zone))
    coord['zone'] =1
    #Pour plusieurs zones modification TH
    # if nb_zone >= 2:
    #     for i in range(2,int(nb_zone)+1):
    #         coord['zone'] = np.where(coord['z'] <= coord.loc[((i-1)*z_top*100/nb_zone),'z'], i,coord['zone'])
            
    
    #coord['zone'] = np.where(coord['z'] <= coord.loc[i-1,'z'], i,coord['zone'])
            #coord['zone'] = np.where(coord['z'] <= thk2, 2,coord['zone'])

    #Write new ginette files
    coord.zone.to_csv(f_zone, index = False, header=False)
    
    # close files    
    f_zone.close()
    f_coor.close()
    
    ##Modifier le fichier de paramètres en fonction du nombre de couches voulues
    f=open('E_zone_parameter_backup.dat','w')
    for i in range(0,nb_zone):
        f.write(str(i+1)+"	[k1]	[n1] [a1] [nVG1] [swres1] [l"+str(i)+"] [c"+str(i)+"] [r1]\n")
    f.close()
    
    return nb_cell

def ginettout(c,l,flux,water_table,nb_zone,alpha,n,s):

    # constant parameters
    ## intrinsic permeability [m2]  k=K*mu/(rho*g)
    ## K hydraulic conductivity [m.s-1]
    ## mu viscosity [Pa.s]
    ## rho density [kg.m-3]
    ## g gravity  9.81 [m2.s-1]
    val_k=4e-15
    # solid grain density rho_s=val_r  [kg.m-3]
    val_r=900

    # porosity
    REF_n=0.10 # \Phi

    # Van Genuchten parameters
    # REF_a=1.5000 #m-1 alpha_vg
    REF_a=alpha
    # REF_nVG= 1.80  # n_vg
    REF_nVG=n
    # REF_swres=0.1 # S_wr
    REF_swres=s

    # Heat capacity is calculated  by the following relationship
    #  c_pm= c_w r_w n  * sat+ c_s r (1-n) + c_a r_a n * (1-sat)
    # density
    # c_s solid specific heat capacity
    #val_c= c_s m2/s2/C I advice to let this value constant.
    # There are no way to calibrate the both parameter rho and c in the same time.
    #c_w=4185D+00	       m2/s2/C
    #r_w=1000  kg/m3    subprocess.call(["./ginette"])

    # solid density r=val_r 
    REF_r=1164

    # REF_c1=c
    # REF_c2=c
    # REF_l=0.7 #

    # Boundary condition water level (meter)
    flux_sup=flux
    REF_WT=water_table
    
    # Initial conditions 
    f_IC_bck=open("E_cdt_initiale_backup.dat","r")
    IC_model=f_IC_bck.read()
    IC_model=IC_model.replace('[head_ini]', '%05.2fD+00' % REF_WT)
    ## write the boundary conditions
    f_bc_bck=open("E_cdt_aux_limites_backup.dat", "r")
    bc_model=f_bc_bck.read()
    bc_model=bc_model.replace('[top]', '%08.2fD+00' % flux_sup)
    bc_model=bc_model.replace('[bot]','%08.2fD+00' % REF_WT)
    # bc_model=bc_model.replace('[dt]','%006.0fD+00' % dt)

    
    ########### Zone of parameters
    f_paramZ_bck=open("E_zone_parameter_backup.dat", "r")
    f_paramZ_new = open("E_zone_parameter.dat", 'w')
    f_bc_new = open("E_cdt_aux_limites.dat", 'w')
    f_IC_new=open("E_cdt_initiale.dat","w")
    param_zone=f_paramZ_bck.read()
    
    # replace the parameter values
    param_zone=param_zone.replace('[k1]','%8.2e' % val_k)
    param_zone=param_zone.replace('[n1]','%6.2f' % REF_n)
    param_zone=param_zone.replace('[r1]','%6.2f' % val_r)
    param_zone=param_zone.replace('[a1]','%8.2f' % REF_a)
    param_zone=param_zone.replace('[nVG1]','%6.2f' % REF_nVG)
    
    param_zone=param_zone.replace('[r1]','%6.2f' % REF_r)
    param_zone=param_zone.replace('[swres1]','%6.2f' % REF_swres)
    
    for i in range(0,nb_zone):
        param_zone=param_zone.replace('[c'+str(i)+']','%6.2f' % c[i])
        param_zone=param_zone.replace('[l'+str(i)+']','%6.2f' % l[i])
        #param_zone=param_zone.replace('[swres'+str(i)+']','%6.2f' % swres[i])

    
    # param_zone=param_zone.replace('[k1]','%8.2e' % val_k)
    # param_zone=param_zone.replace('[n1]','%6.2f' % REF_n)
    # param_zone=param_zone.replace('[r1]','%6.2f' % val_r)
    # param_zone=param_zone.replace('[a1]','%8.2f' % REF_a)
    # param_zone=param_zone.replace('[nVG1]','%6.2f' % REF_nVG)
    # param_zone=param_zone.replace('[swres1]','%6.2f' % REF_swres)
    # param_zone=param_zone.replace('[l1]','%6.2f' % REF_l)
    # param_zone=param_zone.replace('[c1]','%6.2f' % REF_c[0])
    # param_zone=param_zone.replace('[r1]','%6.2f' % REF_r)
    
    # param_zone=param_zone.replace('[k2]','%8.2e' % val_k)
    # param_zone=param_zone.replace('[n2]','%6.2f' % REF_n)
    # param_zone=param_zone.replace('[r2]','%6.2f' % val_r)
    # param_zone=param_zone.replace('[a2]','%8.2f' % REF_a)
    # param_zone=param_zone.replace('[nVG2]','%6.2f' % REF_nVG)
    # param_zone=param_zone.replace('[swres2]','%6.2f' % REF_swres)
    # param_zone=param_zone.replace('[l2]','%6.2f' % REF_l)
    # param_zone=param_zone.replace('[c2]','%6.2f' % REF_c2)
    # param_zone=param_zone.replace('[r2]','%6.2f' % REF_r)
    
    
    #Write new ginette files
    f_IC_new.write(IC_model)
    f_paramZ_new.write(param_zone)
    f_bc_new.write(bc_model)
        
    # close files    
    f_paramZ_new.close()
    f_paramZ_bck.close()

    f_IC_new.close()
    f_bc_bck.close()
    f_bc_new.close()


    # run ginette
    return os.system("timeout 6 ./ginette")
    # subprocess.call(["./ginette"])