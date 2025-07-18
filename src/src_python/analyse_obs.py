import pandas as pd
import numpy as np
import os
from glob import glob

# === FONCTION INDICATEUR DÉJÀ DÉFINIE === (cf. messages précédents)
def calcul_indicateurs(...):  # comme avant
    pass

# === DÉTECTION ÉVÉNEMENTS PLUIE (simple) ===
def detecter_evenements(df, var='pluie', seuil=2):
    if var not in df.columns:
        return []
    return df[df[var] > seuil].index

# === TRAITEMENT MULTI-FICHIERS ===
dossier = './donnees/'
fichiers = glob(os.path.join(dossier, '*.csv'))

resultats = []

indicateurs = []

for t in crues.index:  # ou événements caniculaires
    res = calcul_indicateurs(
        df_smooth, t_event=t,
        var_forcage=col_pluie,
        var_reponse=col_hauteur_eau,  # ou col_temp_eau
        window='2D'
    )
    if res:
        indicateurs.append(res)

df_indicateurs = pd.DataFrame(indicateurs)


for fichier in fichiers:
    # === Lecture du fichier ===
    df = pd.read_csv(fichier, parse_dates=['datetime'])
    df.set_index('datetime', inplace=True)
    df = df.resample('15min').mean().interpolate('time')
    
    # === Identifier le type ===
    nom = os.path.basename(fichier).replace('.csv','')
    if 'PIEZO' in nom:
        type_station = 'piezometre'
        var_reponse = 'temp_eau' if 'temp' in df.columns else 'niveau'
    elif '_H' in nom:
        type_station = 'hauteur_riviere'
        var_reponse = 'hauteur'
    elif '_Q' in nom:
        type_station = 'debit_riviere'
        var_reponse = 'debit'
    else:
        continue  # inconnu

    # === Détection événements ===
    t_events = detecter_evenements(df, var='pluie', seuil=2)

    for t in t_events:
        indic = calcul_indicateurs(
            df, t_event=t,
            var_forcage='pluie',
            var_reponse=var_reponse,
            window='2D'
        )
        if indic:
            indic['station'] = nom
            indic['type'] = type_station
            resultats.append(indic)

# === Résultat final ===
df_resultats = pd.DataFrame(resultats)
df_resultats.to_csv("indicateurs_reactivite.csv", index=False)
