# -*- coding: utf-8 -*-
"""
Outils communs aux cas de test InterFrost 2D de Ginette (TH2 "inclusion gelée"
et TH3 "talik", Grenier et al. 2018, Adv. Water Resour. 114, 196-218) :

- lecture des mesures de performance écrites par Ginette (ytest=TH2/TH3) dans
  S_TH2.dat / S_TH3.dat (Fortran unformatted sequential, un enregistrement de
  real*4 par pas de sortie) et du champ final S_pts ;
- conversion des résultats de référence InterFrost (classeurs Excel du site
  wiki.lsce.ipsl.fr/interfrost, un onglet par mesure et par gradient, 13
  modèles anonymes) en CSV long compact, et lecture de ces CSV ;
- comparaison d'une courbe Ginette à l'enveloppe des modèles de référence.
"""
import numpy as np
import pandas as pd

# Gradients de charge de l'intercomparaison (fiches TH2/TH3 + complément du
# 11/05/2015) ; l'indice GrH1..4 des onglets Excel suit cet ordre.
GRADIENTS = {"TH2": [0.0, 0.03, 0.09, 0.15], "TH3": [0.03, 0.06, 0.09, 0.15]}

# Colonnes des enregistrements binaires de Ginette (après le temps)
PM_COLUMNS = {
    "TH2": ["t_s", "PM1_Tmin_C", "PM3_Vwater_m3", "PM2_Jout_W"],
    "TH3": ["t_s", "PM1_Keq_m_s", "PM2_Jcond_W", "PM3_Jsens_J", "PM4_Pt1_C", "PM4_Pt2_C"],
}

# Correspondance onglet Excel -> nom de mesure (le classeur TH2 écrit "MP2"/"MP3")
SHEETS = {
    "TH2": {"PM1": "TH2_PM1_GrH{k}", "PM2": "TH2_MP2_GrH{k}", "PM3": "TH2_MP3_GrH{k}"},
    # TH3_PM2 (flux conductif aux frontières haut/bas) a été abandonné par
    # l'intercomparaison ("large sensitivity to spatial discretization") : pas
    # de courbes de référence.
    "TH3": {"PM1": "TH3_PM1_GrH{k}", "PM3": "TH3_MP3_GrH{k}",
            "PM4_Pt1": "TH3_MP4_Pt1_GrH{k}", "PM4_Pt2": "TH3_MP4_Pt2_GrH{k}"},
}


def read_pm_records(path, case):
    """Lit S_TH2.dat / S_TH3.dat : enregistrements Fortran unformatted
    sequential de n real*4 (4 octets de marqueur avant et après chacun).
    Renvoie un DataFrame avec les colonnes PM_COLUMNS[case]."""
    cols = PM_COLUMNS[case]
    n = len(cols)
    raw = np.fromfile(path, dtype=np.uint8)
    rec = 8 + 4 * n
    nrec = len(raw) // rec
    body = raw[:nrec * rec].reshape(nrec, rec)[:, 4:4 + 4 * n].copy()
    vals = np.frombuffer(body.tobytes(), dtype=np.float32).reshape(nrec, n).astype(float)
    return pd.DataFrame(vals, columns=cols)


def read_s_pts_2d(path, nm):
    """Champ final écrit dans S_pts (ytest=TH2/TH3) : 1 enregistrement real(t),
    puis nm enregistrements (real(pr), real(temp), real(sw)). Renvoie
    (t, pr, temp, sw), tableaux de longueur nm dans l'ordre des mailles."""
    raw = np.fromfile(path, dtype=np.uint8)
    assert len(raw) == 12 + nm * 20, f"S_pts : {len(raw)} octets, attendu {12 + nm*20}"
    t = float(np.frombuffer(raw[4:8].tobytes(), dtype=np.float32)[0])
    body = np.frombuffer(raw[12:].tobytes(), dtype=np.float32).reshape(nm, 5)
    return t, body[:, 1].astype(float), body[:, 2].astype(float), body[:, 3].astype(float)


def read_coordinates(path, nm):
    """E_coordonnee.dat écrit par Ginette (maillage automatique) : i, x, z."""
    c = pd.read_csv(path, sep=r"\s+", header=None, names=["i", "x", "z"]).tail(nm)
    return c["x"].to_numpy(), c["z"].to_numpy()


def convert_reference_xlsx(xlsx, case, out_csv, max_points=300):
    """Convertit le classeur InterFrost (resultsinterfrost_th2.xlsx / _th3.xlsx)
    en CSV long : case, pm, gradient, model, t_s, value. Chaque courbe est
    sous-échantillonnée à max_points points (uniformément en indice) pour
    garder le CSV léger ; les courbes ont de 25 à 20000 points."""
    x = pd.ExcelFile(xlsx)
    rows = []
    for pm, pattern in SHEETS[case].items():
        for k, grad in enumerate(GRADIENTS[case], start=1):
            sheet = pattern.format(k=k)
            if sheet not in x.sheet_names:
                continue
            d = x.parse(sheet)
            for m in range(1, 40):
                tc, pc = f"Model {m} - Time", f"Model {m} - PM"
                if tc not in d.columns or pc not in d.columns:
                    continue
                t = pd.to_numeric(d[tc], errors="coerce")
                v = pd.to_numeric(d[pc], errors="coerce")
                ok = t.notna() & v.notna()
                if ok.sum() < 2:
                    continue
                t, v = t[ok].to_numpy(), v[ok].to_numpy()
                order = np.argsort(t, kind="stable")
                t, v = t[order], v[order]
                if len(t) > max_points:
                    idx = np.unique(np.round(np.linspace(0, len(t) - 1, max_points)).astype(int))
                    t, v = t[idx], v[idx]
                rows.append(pd.DataFrame({"case": case, "pm": pm, "gradient": grad,
                                          "model": m, "t_s": t, "value": v}))
    ref = pd.concat(rows, ignore_index=True)
    ref.to_csv(out_csv, index=False, float_format="%.6g")
    return ref


def load_reference(csv, pm, gradient):
    """Courbes de référence d'une mesure et d'un gradient : {model: (t, v)}."""
    ref = pd.read_csv(csv)
    sel = ref[(ref["pm"] == pm) & (np.isclose(ref["gradient"], gradient))]
    return {m: (g["t_s"].to_numpy(), g["value"].to_numpy()) for m, g in sel.groupby("model")}


def envelope(curves, t):
    """Min, médiane et max des courbes de référence interpolées aux temps t
    (chaque courbe n'est utilisée que sur son propre intervalle de temps)."""
    stack = []
    for tm, vm in curves.values():
        v = np.interp(t, tm, vm, left=np.nan, right=np.nan)
        stack.append(v)
    stack = np.array(stack)
    with np.errstate(all="ignore"):
        return (np.nanmin(stack, axis=0), np.nanmedian(stack, axis=0), np.nanmax(stack, axis=0),
                np.sum(~np.isnan(stack), axis=0))


def threshold_time(t, v, level=0.0, rising=True):
    """Premier instant où v franchit level (interpolation linéaire), NaN sinon."""
    t, v = np.asarray(t, float), np.asarray(v, float)
    s = (v - level) if rising else (level - v)
    idx = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    return t[i] + (t[i + 1] - t[i]) * (-s[i]) / (s[i + 1] - s[i])


def compare_to_reference(t, v, curves):
    """Écart d'une courbe Ginette à l'enveloppe de référence : RMSE à la
    médiane, fraction du temps passée dans [min, max], et écart normalisé par
    la largeur de l'enveloppe (moyenne de |v - médiane| / (max - min))."""
    lo, med, hi, n = envelope(curves, t)
    ok = n >= 3
    if ok.sum() == 0:
        return dict(rmse_median=np.nan, frac_in_envelope=np.nan, rel_spread=np.nan)
    inside = (v[ok] >= lo[ok] - 1e-12) & (v[ok] <= hi[ok] + 1e-12)
    width = np.maximum(hi[ok] - lo[ok], 1e-12)
    return dict(rmse_median=float(np.sqrt(np.mean((v[ok] - med[ok]) ** 2))),
                frac_in_envelope=float(inside.mean()),
                rel_spread=float(np.mean(np.abs(v[ok] - med[ok]) / width)))
