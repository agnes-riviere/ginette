
def one_set(dt,nb_day,z_bottom,dz,nb_cell,cell1,cell2,cell3,cell4,val_k,val_n,val_l,val_r):
    import os
    import numpy as np
    from pathlib import Path
    import pandas as pd
    from scipy import interpolate
    import matplotlib.pyplot as plt
    from IPython.display import display
    import subprocess
    import shutil
    state=0

    setup_model=setup_model.replace('[dt]','%06.0fD+00' % dt)
    setup_model=setup_model.replace('[nb_day]','%06.0f' % nb_day)
    setup_model=setup_model.replace('[z_bottom]','%6.2e' % z_bottom)
    setup_model=setup_model.replace('[az]','%7.3e' % -z_bottom)
    setup_model=setup_model.replace('[dz]','%6.2e' % dz)
    setup_model=setup_model.replace('[nb_cell]','%05.0f' % nb_cell)
    setup_model=setup_model.replace('[cell1]','%05d' % cell1)
    setup_model=setup_model.replace('[cell2]','%05d' % cell2)
    setup_model=setup_model.replace('[cell3]','%05d' % cell3)
    setup_model=setup_model.replace('[cell4]','%05d' % cell4)

    setup_model=setup_model.replace('[state]','%1i' % state)

    ########### Zone
    f_coor=open("E_coordonnee.dat", "r")
    f_zone=open("E_zone.dat", 'w')
    f_paramZ_bck=open("E_zone_parameter_backup.dat", "r")
    f_cdi_bck = open("E_cdt_initiale_bck.dat", 'r')
    f_cdl_bck = open("E_cdt_aux_limites_bck.dat",'r')
    f_paramZ_new = open("E_zone_parameter.dat", 'w')
    f_param_new = open("E_parametre.dat", 'w')
    f_cdi_new = open("E_cdt_initiale.dat", 'w')
    f_cdl_new = open("E_cdt_aux_limites.dat",'w')
    param_zone=f_paramZ_bck.read()
    cdi=f_cdi_bck.read()
    cdl=f_cdl_bck.read()
    coord=pd.DataFrame()
    coord = pd.read_csv(f_coor, names=["id", "x", "z"], header=None, delim_whitespace=True)



    param_zone=param_zone.replace('[k1]','%8.2e' % val_k)
    param_zone=param_zone.replace('[n1]','%6.2f' % val_n)
    param_zone=param_zone.replace('[l1]','%6.2f' % val_l)
    param_zone=param_zone.replace('[r1]','%6.2f' % val_r)


    coord['zone'] =1

    if nb_zone >= 2:
        coord['zone'] = np.where(coord['z'] <= thk2, 2,coord['zone'])
        coord['zone'] = np.where(coord['z'] <= thk3, 3,coord['zone'])

    #display(coord)

    f_param_new.write(setup_model)
    f_paramZ_new.write(param_zone)
    coord.zone.to_csv(f_zone, index = False, header=False)
    f_zone.close()
    f_paramZ_new.close()
    f_paramZ_bck.close()

    f_coor.close()

    print("Ginette's steady")

    # change parameters in input files for steady state
    setup_model=setup_model.replace('[state]','%1i' % state)
    ichi2=0
    itempi=0
    iclchgt=0
    cdi=cdi.replace('[ichi2]','%1i' % ichi2)
    cdi=cdi.replace('[itempi]','%1i' % itempi)
    cdl=cdl.replace('[iclchgt]','%1i' % iclchgt)
    f_cdl_new.write(cdl)
    f_cdi_new.write(cdi)

    f_param_bck.close()
    f_param_new.close()
    f_cdi_bck.close()
    f_cdi_new.close()
    f_cdl_bck.close()
    f_cdl_new.close()

    subprocess.call(["./ginette"])   


    print("Ginette transient")
    f_param_bck2=open("E_parametre_backup.dat", "r")
    f_cdi_bck2 = open("E_cdt_initiale_bck.dat", 'r')
    f_cdl_bck2 = open("E_cdt_aux_limites_bck.dat",'r')
    f_param_new2 = open("E_parametre.dat", 'w')
    f_cdi_new2 = open("E_cdt_initiale.dat", 'w')
    f_cdl_new2 = open("E_cdt_aux_limites.dat",'w')
    cdi2=f_cdi_bck2.read()
    cdl2=f_cdl_bck2.read()
    setup_model2=f_param_bck2.read()
    setup_model2=setup_model2.replace('[dt]','%06.0fD+00' % dt)
    setup_model2=setup_model2.replace('[nb_day]','%06.0f' % nb_day)
    setup_model2=setup_model2.replace('[z_bottom]','%6.2e' % z_bottom)
    setup_model2=setup_model2.replace('[az]','%7.3e' % -z_bottom)
    setup_model2=setup_model2.replace('[dz]','%6.2e' % dz)
    setup_model2=setup_model2.replace('[nb_cell]','%05.0f' % nb_cell)
    setup_model2=setup_model2.replace('[cell1]','%05d' % cell1)
    setup_model2=setup_model2.replace('[cell2]','%05d' % cell2)
    setup_model2=setup_model2.replace('[cell3]','%05d' % cell3)
    setup_model2=setup_model2.replace('[cell4]','%05d' % cell4)

    state=1
    ichi2=1
    itempi=1
    iclchgt=1
    cdi2=cdi2.replace('[ichi2]','%1i' % ichi2)
    cdi2=cdi2.replace('[itempi]','%1i' % itempi)
    cdl2=cdl2.replace('[iclchgt]','%1i' % iclchgt)
    setup_model2=setup_model2.replace('[state]','%1i' % state)

    f_cdl_new2.write(cdl2)
    f_cdi_new2.write(cdi2)
    f_param_new2.write(setup_model2)

    # paste pressure results from steady state
    shutil.copyfile('S_pression_charge_temperature.dat','E_pression_initiale.dat')

    f_param_bck2.close()
    f_param_new2.close()
    f_cdi_bck2.close()
    f_cdi_new2.close()
    f_cdl_bck2.close()
    f_cdl_new2.close()

    subprocess.call(["./ginette"])   

    os.chdir(path_mini_lomos)# use copyfile()
