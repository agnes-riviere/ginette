import os

from pathlib import Path
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib as mpl
from IPython.display import display
import subprocess
import glob
import itertools
import seaborn as sns
import array

def dic_setup_model(params=None):
    # Définir les valeurs par défaut
    default_params = {
        "dt": 900,
        "nb_day": 30*7,
        "date_simul_bg": "2022-05-01 14:00:00",
        "state": 1,
        "z_top": -0.1,
        "z_bottom": -0.4,
        "dz": 0.01,
        "nb_zone": 2
    }
# Si des paramètres sont fournis, mettre à jour les valeurs par défaut
    if params is not None:
        for key in default_params.keys():
            if key not in params:
                params[key] = default_params[key]

    return params

def setup_model(dt=None, nb_day=None, date_simul_bg=None, state=None, 
                z_top=None, z_bottom=None, dz=None, nb_zone=None):
    # Créer un dictionnaire avec les paramètres fournis
    params = {"dt": dt, "nb_day": nb_day, "date_simul_bg": date_simul_bg, "state": state, 
              "z_top": z_top, "z_bottom": z_bottom, "dz": dz, "nb_zone": nb_zone}

    # Obtenir les paramètres par défaut pour ceux qui n'ont pas été fournis
    params = dic_setup_model(params)
    date_simul_bg = pd.to_datetime(date_simul_bg)
    az = abs(z_top-z_bottom)
    nb_cell = az/dz



    print("la simulation commence à", date_simul_bg)

    with open("E_parametre_backup.dat", "r") as f_param_bck, open("E_parametre.dat", 'w') as f_param_new:
        setup_model = f_param_bck.read()
        setup_model = setup_model.replace('[dt]', '%06.0fD+00' % dt)
        setup_model = setup_model.replace('[state]', '%1i' % state)
        setup_model = setup_model.replace('[nb_day]', '%06.0f' % nb_day)
        setup_model = setup_model.replace('[z_top]', '%7.2e' % z_top)
        setup_model = setup_model.replace('[z_bottom]', '%7.2e' % z_bottom)
        setup_model = setup_model.replace('[az]', '%7.3e' % az)
        setup_model = setup_model.replace('[dz]', '%6.2e' % dz)
        setup_model = setup_model.replace('[nb_cell]', '%05.0f' % nb_cell)

        Obs1 = z_top
        Obs2 = z_top-0.1
        Obs3 = z_top-0.2
        Obs4 = z_top-0.3

        cell1 = abs((z_top-Obs1)/dz)+1
        cell2 = abs((z_top-Obs2)/dz)
        cell3 = abs((z_top-Obs3)/dz)
        cell4 = abs((z_top-Obs4)/dz)

        setup_model = setup_model.replace('[cell1]', '%05d' % cell1)
        setup_model = setup_model.replace('[cell2]', '%05d' % cell2)
        setup_model = setup_model.replace('[cell3]', '%05d' % cell3)
        setup_model = setup_model.replace('[cell4]', '%05d' % cell4)

        f_param_new.write(setup_model)