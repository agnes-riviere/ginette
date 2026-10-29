# -*- coding: utf-8 -*-
"""
Estimation de la diffusivité thermique et de la vitesse de Darcy à partir de profils
verticaux de température, par la méthode de Stallman (1965) / Tabbagh et al.,
reprise notamment dans Cheviron, Guérin, Tabbagh, Bendjoudi (2005),
"Determining long-term effective groundwater recharge by analyzing vertical soil
temperature profiles at meteorological stations", Water Resources Research 41, W09501.

Principe :
----------
Pour un milieu homogène soumis à une variation de température sinusoïdale en
surface (période P, pulsation w=2*pi/P) et un écoulement vertical de Darcy uz,
la solution de Stallman (1965) de l'équation de la chaleur avec advection est :

    T(z, t) = T0 + A0 * exp(-a*z) * cos(w*t - b*z)

où a (amortissement, rad/m... en fait 1/m) et b (déphasage, rad/m) dépendent de
la diffusivité thermique Gamma = k/Cv [m2/s] et de la vitesse frontale thermique
v = uz * Cw / Cv [m/s] (Cw = chaleur volumique de l'eau, Cv = chaleur volumique
du milieu saturé) :

    lambda = a + i*b = (-v + sqrt(v**2 + 4*i*Gamma*w)) / (2*Gamma)

En l'absence d'écoulement (v=0), a=b=sqrt(w/(2*Gamma)) : c'est la base des
"diffusivités apparentes" Gamma_ampl et Gamma_phase couramment calculées
indépendamment à partir de l'amortissement d'amplitude et du déphasage entre
deux profondeurs.

En présence d'écoulement, a et b divergent et on peut inverser le système pour
obtenir Gamma et v simultanément (voir invert_v_gamma ci-dessous, dérivé de la
relation de Stallman ci-dessus et validé numériquement par un aller-retour
direct/inverse sur données synthétiques).
"""
import numpy as np
import pandas as pd


def fit_sinusoid(t_seconds, values, period_seconds, trend_order=3):
    """
    Ajuste par moindres carrés un modèle T(t) = trend(t) + A*cos(w*t) + B*sin(w*t),
    où trend(t) est un polynôme de degré `trend_order` (par défaut cubique, comme
    dans la note de Tabbagh pour absorber les variations lentes/dérives).

    Parameters:
    - t_seconds: temps en secondes depuis un instant de référence commun (array-like).
    - values: valeurs de température correspondantes (array-like), NaN autorisés (ignorés).
    - period_seconds: période de l'oscillation à extraire (86400 pour le cycle diurne).
    - trend_order: degré du polynôme de tendance (0 = juste une moyenne).

    Returns:
    - dict avec 'amplitude' (>=0), 'phase' (rad, tel que le signal ~ amplitude*cos(w*t-phase)),
      'mean', et 'residual_std' (écart-type des résidus, indicateur de qualité de l'ajustement).
    """
    t = np.asarray(t_seconds, dtype=float)
    y = np.asarray(values, dtype=float)
    mask = np.isfinite(t) & np.isfinite(y)
    t, y = t[mask], y[mask]

    w = 2 * np.pi / period_seconds
    # Le polynôme de tendance est évalué sur un temps recentré/mis à l'échelle (en jours)
    # pour garder la matrice de régression bien conditionnée (t**3 en secondes exploserait
    # numériquement sur une fenêtre de plusieurs dizaines de jours). La phase et l'amplitude
    # du terme harmonique, elles, utilisent bien w*t en secondes (temps physique réel).
    t_trend = (t - t.mean()) / 86400.0
    cols = [t_trend**k for k in range(trend_order + 1)]
    cols += [np.cos(w * t), np.sin(w * t)]
    X = np.column_stack(cols)

    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    A, B = coeffs[-2], coeffs[-1]
    amplitude = np.hypot(A, B)
    phase = np.arctan2(B, A)  # tel que A*cos(wt)+B*sin(wt) = amplitude*cos(wt-phase)

    residuals = y - X @ coeffs
    return {
        'amplitude': amplitude,
        'phase': phase,
        'mean': coeffs[0],
        'residual_std': float(np.std(residuals)),
    }


def apparent_diffusivities(amp_shallow, amp_deep, phase_shallow, phase_deep, dz, period_seconds):
    """
    Calcule les diffusivités thermiques apparentes Gamma_ampl et Gamma_phase entre
    deux profondeurs, en supposant un milieu homogène SANS écoulement (v=0).
    Ce sont les grandeurs tabulées dans la note de Tabbagh (colonnes Gamma_ampl / Gamma_phase).

    Parameters:
    - amp_shallow, amp_deep: amplitudes du signal harmonique aux deux profondeurs (deep = plus profonde).
    - phase_shallow, phase_deep: phases (rad) aux deux profondeurs (voir fit_sinusoid).
    - dz: écart de profondeur (m), toujours positif.
    - period_seconds: période de l'oscillation utilisée (86400 pour le cycle diurne).

    Returns:
    - (Gamma_ampl, Gamma_phase) en m2/s. Le déphasage est "unwrappé" en supposant
      qu'il croît de façon monotone avec la profondeur (retard croissant).
    """
    w = 2 * np.pi / period_seconds
    a = np.log(amp_shallow / amp_deep) / dz
    # Déphasage : on force un retard positif (le signal profond est en retard sur le signal de surface)
    dphi = (phase_shallow - phase_deep) % (2 * np.pi)
    b = dphi / dz

    gamma_ampl = w * dz**2 / (2 * (a * dz)**2) if a != 0 else np.nan
    gamma_phase = w / (2 * b**2) if b != 0 else np.nan
    return gamma_ampl, gamma_phase


def invert_v_gamma(a, b, period_seconds):
    """
    Inversion jointe de la diffusivité thermique Gamma et de la vitesse frontale
    thermique v (= uz*Cw/Cv) à partir de l'amortissement a=ln(Ampl_haut/Ampl_bas)/dz
    et du déphasage b=Dphi/dz mesurés entre deux profondeurs, en résolvant
    exactement la relation de Stallman (voir docstring du module).

    Validé numériquement (aller-retour exact sur données synthétiques, y compris
    v négatif = exfiltration et v=0 = cas purement diffusif où a=b).

    Parameters:
    - a: taux d'amortissement de l'amplitude (1/m), toujours positif.
    - b: taux de déphasage (rad/m), toujours positif.
    - period_seconds: période de l'oscillation utilisée.

    Returns:
    - (v, Gamma) : v en m/s (vitesse frontale thermique), Gamma en m2/s.
    """
    w = 2 * np.pi / period_seconds
    p = w / b
    gamma = a * p**3 / (w**2 + a**2 * p**2)
    v = p - 2 * gamma * a
    return v, gamma


def darcy_velocity_from_v(v, Cv, Cw=4.18e6):
    """
    Convertit la vitesse frontale thermique v (m/s) en vitesse de Darcy uz (m/s),
    via v = uz * Cw / Cv (Cw = chaleur volumique de l'eau, Cv = chaleur volumique
    du milieu saturé).
    """
    return v * Cv / Cw


def solid_thermal_conductivity(k_bulk, porosity, kw=0.598):
    """
    Convertit la conductivité thermique du milieu SATURÉ (bulk, k_bulk = Gamma*Cv,
    ce que la méthode de Stallman/Cheviron donne directement) en conductivité
    thermique de la fraction SOLIDE seule (ks) — c'est ks, pas k_bulk, qui
    correspond au paramètre `lam` calibré dans le grid search Ginette
    (1_define_grid_search.py) : lam ne dépend que du milieu poreux (composition
    minérale), pas de sa saturation en eau.

    IMPORTANT : la loi de mélange utilisée ici DOIT être la même que celle que
    Ginette applique en interne pour recombiner ks et kw en k_bulk, sinon les
    bornes de lam dérivées ne correspondent pas à ce que le modèle fera
    réellement de ces valeurs. Ginette sélectionne cette loi via le mot-clé
    `ymoycondtherm` de E_p_therm.dat (voir src/ginette_V2.f90 ~ligne 1907) :
    ARITH (moyenne arithmétique, réglage par défaut), GEOME (moyenne
    géométrique/produit, utilisée par Cheviron et al. 2005 éq. 8 et la note de
    Tabbagh), WOODS, LUNAR, NEUMA. Avec ymoycondtherm=ARITH :

        k_bulk = kw*n + ks*(1-n)    =>    ks = (k_bulk - kw*n) / (1 - n)

    Si E_p_therm.dat est reconfiguré avec un autre ymoycondtherm, cette
    fonction doit être mise à jour en conséquence (ex: modèle produit ci-dessus
    en commentaire, utilisé par une version précédente de ce module).

    Parameters:
    - k_bulk: conductivité thermique du milieu saturé [W/m/K] (= Gamma * Cv).
    - porosity: porosité n (0 à 1).
    - kw: conductivité thermique de l'eau [W/m/K] (0.598, alandae dans Ginette).
    """
    return (k_bulk - kw * porosity) / (1 - porosity)


def water_viscosity(temperature_celsius):
    """
    Viscosité dynamique de l'eau [Pa.s] en fonction de la température [°C],
    via la relation de Vogel (précision suffisante entre 0 et 40°C environ,
    largement utilisée en hydrogéologie pour ce type de calcul) :

        mu(T) = 2.414e-5 * 10**(247.8 / (T + 133.15))
    """
    return 2.414e-5 * 10 ** (247.8 / (temperature_celsius + 133.15))


def darcy_law_K(uz, dh, dz):
    """
    Conductivité hydraulique K [m/s] déduite de la loi de Darcy à partir d'une
    vitesse de Darcy uz [m/s] et d'un gradient hydraulique dh/dz mesuré
    indépendamment (ici via une différence de charge rivière-aquifère) :

        uz = K * dh/dz   =>   K = uz / (dh/dz) = uz * dz / dh

    Convention de signe : dh = h_amont - h_aval (ex: h_rivière - h_piézomètre),
    dz = distance entre les deux points de mesure (toujours positif). dh>0 (la
    rivière est plus haute que l'aquifère) doit correspondre à une infiltration
    (uz>0 avec la convention "z positif vers le bas" utilisée dans ce module).
    """
    return uz * dz / dh


def intrinsic_permeability(K, temperature_celsius, rho=1000.0, g=9.81):
    """
    Convertit une conductivité hydraulique K [m/s] en perméabilité intrinsèque
    k [m2] (c'est k, pas K, qui est calibré dans le grid search Ginette, cf.
    log_k = log10(k) dans 1_define_grid_search.py) :

        K = k * rho * g / mu(T)   =>   k = K * mu(T) / (rho * g)

    La viscosité de l'eau dépend de la température : on utilise ici la
    température de l'eau du ru mesurée (temperature_eau / TempMolo), moyennée
    sur la même fenêtre que le reste de l'analyse.
    """
    mu = water_viscosity(temperature_celsius)
    return K * mu / (rho * g)


def analyze_depth_pairs(fits, depths, period_seconds, Cv, Cw=4.18e6):
    """
    À partir des ajustements sinusoïdaux (sortie de fit_sinusoid) pour un ensemble
    de profondeurs, calcule pour chaque paire (profondeur peu profonde, profondeur
    plus profonde) : Gamma_ampl, Gamma_phase, la vitesse de Darcy inversée (uz) et
    un indicateur de cohérence (accord entre Gamma_ampl et Gamma_phase, cf. note
    de Tabbagh : un fort désaccord signale un milieu non homogène ou une advection
    importante à interpréter avec prudence).

    Parameters:
    - fits: dict {nom_profondeur: résultat de fit_sinusoid(...)}.
    - depths: dict {nom_profondeur: profondeur en m (valeur positive = sous la surface)}.
    - period_seconds: période utilisée pour les ajustements (doit être cohérente avec fits).
    - Cv: chaleur volumique du milieu saturé [J/m3/K] (à adapter selon porosité/lithologie).
    - Cw: chaleur volumique de l'eau [J/m3/K] (4.18e6 par défaut).

    Returns:
    - pd.DataFrame avec une ligne par paire de profondeurs.
    """
    names = sorted(depths, key=lambda n: depths[n])  # du plus proche de la surface au plus profond
    rows = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            dz = depths[n2] - depths[n1]
            f1, f2 = fits[n1], fits[n2]

            gamma_ampl, gamma_phase = apparent_diffusivities(
                f1['amplitude'], f2['amplitude'], f1['phase'], f2['phase'], dz, period_seconds
            )

            a = np.log(f1['amplitude'] / f2['amplitude']) / dz
            dphi = (f1['phase'] - f2['phase']) % (2 * np.pi)
            b = dphi / dz
            v, gamma_joint = invert_v_gamma(a, b, period_seconds)
            uz = darcy_velocity_from_v(v, Cv, Cw)

            coherence = abs(gamma_ampl - gamma_phase) / gamma_phase if gamma_phase else np.nan

            rows.append({
                'pair': f'{n1}-{n2}',
                'dz_m': dz,
                'gamma_ampl_1e-6': gamma_ampl * 1e6,
                'gamma_phase_1e-6': gamma_phase * 1e6,
                'coherence_rel_diff': coherence,
                'gamma_joint_1e-6': gamma_joint * 1e6,
                'v_thermal_m_s': v,
                'uz_darcy_mm_day': uz * 1000 * 86400,
                'uz_darcy_mm_yr': uz * 1000 * 86400 * 365.25,
            })
    return pd.DataFrame(rows)


def decompose_signal(t_seconds, values, period_seconds, trend_order=3):
    """
    Comme fit_sinusoid, mais renvoie en plus la décomposition temporelle complète
    (tendance, harmonique diurne, résidu large-bande) point par point, nécessaire
    pour reconstruire un signal corrigé (voir predict_true_surface_signal)
    plutôt que juste amplitude/phase globales sur la fenêtre.
    """
    t = np.asarray(t_seconds, dtype=float)
    y = np.asarray(values, dtype=float)

    w = 2 * np.pi / period_seconds
    t_trend = (t - t.mean()) / 86400.0
    n_trend_cols = trend_order + 1
    cols = [t_trend**k for k in range(n_trend_cols)]
    cols += [np.cos(w * t), np.sin(w * t)]
    X = np.column_stack(cols)

    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    A, B = coeffs[-2], coeffs[-1]
    amplitude = np.hypot(A, B)
    phase = np.arctan2(B, A)

    trend = X[:, :n_trend_cols] @ coeffs[:n_trend_cols]
    harmonic = A * np.cos(w * t) + B * np.sin(w * t)
    residual = y - trend - harmonic

    return {
        'amplitude': amplitude, 'phase': phase, 'w': w,
        'trend': trend, 'harmonic': harmonic, 'residual': residual,
    }


def fit_local_attenuation_phase(fits, depths, period_seconds):
    """
    Calcule le taux d'atténuation a(z) [1/m] et de déphasage b(z) [rad/m] LOCAUX
    (pas cumulés depuis la surface) pour chaque paire de profondeurs
    CONSÉCUTIVES, assignés au point milieu de la paire.

    Objectif : vérifier si a,b varient régulièrement (donc de façon extrapolable,
    cf. extrapolate_attenuation_phase) avec la profondeur, ou si une paire donnée
    rompt nettement la tendance des autres — typiquement le signe qu'un des deux
    capteurs de cette paire n'est pas dans le même régime de transport que les
    autres (ex: capteur d'eau libre au-dessus du lit, mélangé par turbulence,
    par opposition à la diffusion/advection dans le sédiment saturé).

    Parameters:
    - fits: dict {nom: résultat de fit_sinusoid(...)}.
    - depths: dict {nom: profondeur en m, positif vers le bas}.
    - period_seconds: période utilisée pour les ajustements.

    Returns:
    - pd.DataFrame avec colonnes ['pair', 'z_mid', 'dz', 'a', 'b'].
    """
    names = sorted(depths, key=lambda n: depths[n])
    rows = []
    for i in range(len(names) - 1):
        n0, n1 = names[i], names[i + 1]
        z0, z1 = depths[n0], depths[n1]
        dz = z1 - z0
        A0, A1 = fits[n0]['amplitude'], fits[n1]['amplitude']
        p0, p1 = fits[n0]['phase'], fits[n1]['phase']
        a = np.log(A0 / A1) / dz
        # Le déphasage attendu entre deux capteurs proches (quelques cm) est petit
        # (bien moins qu'un demi-cycle) : on prend donc l'équivalent de plus petite
        # magnitude de (p1-p0) modulo 2*pi (dans (-pi,pi]), PAS (p0-p1)%(2*pi) qui
        # peut retomber sur le mauvais "tour" (~2*pi moins le vrai petit déphasage)
        # dès que p1>p0, ce qui est le cas normal ici (la phase croît avec la
        # profondeur, cf. depth_profile_regression).
        dphi = (p1 - p0 + np.pi) % (2 * np.pi) - np.pi
        b = dphi / dz
        rows.append({'pair': f'{n0}-{n1}', 'z_mid': (z0 + z1) / 2, 'dz': dz, 'a': a, 'b': b})
    return pd.DataFrame(rows)


def extrapolate_attenuation_phase(local_df, z_target, exclude_pairs=None):
    """
    Régression linéaire de a(z) et b(z) (sortie de fit_local_attenuation_phase)
    sur les paires jugées fiables, puis extrapolation à la profondeur z_target [m].

    Parameters:
    - local_df: sortie de fit_local_attenuation_phase.
    - z_target: profondeur (m) à laquelle extrapoler a et b.
    - exclude_pairs: liste de noms de paires (ex: ['TempMolo-Temp1']) à exclure
      de la régression, typiquement une paire connue pour rompre la tendance.

    Returns:
    - (a_pred, b_pred) : taux d'atténuation [1/m] et de déphasage [rad/m] prédits
      à z_target par continuation de la tendance des paires retenues.
    """
    df = local_df if not exclude_pairs else local_df[~local_df['pair'].isin(exclude_pairs)]
    a_slope, a_intercept = np.polyfit(df['z_mid'], df['a'], 1)
    b_slope, b_intercept = np.polyfit(df['z_mid'], df['b'], 1)
    a_pred = a_slope * z_target + a_intercept
    b_pred = b_slope * z_target + b_intercept
    return a_pred, b_pred


def predict_true_surface_signal(t_seconds, reference_values, anchor_fit, period_seconds,
                                 dz_to_surface, a_pred, b_pred, trend_order=3):
    """
    Corrige un capteur de référence mesuré près de la surface mais PAS au
    contact du sédiment (ex: TempMolo, 1-2cm au-dessus du lit, dans l'eau libre
    du ru) pour estimer la température réelle à la surface du sédiment (z=0).

    Principe : seule la composante diurne (haute fréquence) du capteur de
    référence est remplacée, par celle qu'on prédit en extrapolant À REBOURS,
    sur dz_to_surface, la loi d'atténuation/déphasage interne du sédiment
    (a_pred, b_pred, voir extrapolate_attenuation_phase) à partir d'un capteur
    fiable plus profond (anchor_fit, typiquement Temp1) :

        A0_pred     = A_anchor * exp(a_pred * dz_to_surface)
        phase0_pred = phase_anchor - b_pred * dz_to_surface
        T_surface(t) = tendance_ref(t) + A0_pred*cos(w*t - phase0_pred) + résidu_ref(t)

    La tendance lente et le résidu large-bande du capteur de référence sont
    conservés tels quels (une variation lente s'atténue très peu sur 1-2cm,
    donc pas de raison de la corriger ; le résidu large-bande n'est pas non
    plus concerné par le désaccord observé, qui porte spécifiquement sur
    l'harmonique diurne).

    Parameters:
    - t_seconds: temps (s), commun à reference_values et à l'ajustement anchor_fit.
    - reference_values: série brute du capteur de référence (ex: TempMolo).
    - anchor_fit: sortie de fit_sinusoid() pour le capteur fiable (ex: Temp1).
    - dz_to_surface: distance (m) entre le capteur fiable et la surface (ex: 0.10).
    - a_pred, b_pred: taux d'atténuation/déphasage extrapolés à la surface.

    Returns:
    - dict avec 'corrected' (array, température corrigée à z=0), 'amplitude_pred',
      'phase_pred', 'amplitude_raw'/'phase_raw' (harmonique brute du capteur de
      référence, pour comparaison), et 'trend'/'residual' (composantes réutilisées
      telles quelles).
    """
    decomp = decompose_signal(t_seconds, reference_values, period_seconds, trend_order)
    w = decomp['w']
    t = np.asarray(t_seconds, dtype=float)

    amplitude_pred = anchor_fit['amplitude'] * np.exp(a_pred * dz_to_surface)
    phase_pred = anchor_fit['phase'] - b_pred * dz_to_surface

    harmonic_pred = amplitude_pred * np.cos(w * t - phase_pred)
    corrected = decomp['trend'] + harmonic_pred + decomp['residual']

    return {
        'corrected': corrected,
        'amplitude_pred': amplitude_pred,
        'phase_pred': phase_pred,
        'amplitude_raw': decomp['amplitude'],
        'phase_raw': decomp['phase'],
        'trend': decomp['trend'],
        'residual': decomp['residual'],
    }


def predict_deeper_signal(t_seconds, reference_values, anchor_fit, period_seconds,
                           dz_deeper, a_pred, b_pred, trend_order=3):
    """
    Extrapole le signal du capteur le plus profond disponible (ex: Temp4) vers
    une profondeur supplémentaire dz_deeper où aucune mesure n'existe (ex: la
    condition limite du bas du domaine, décalée de -10cm de plus après
    USE_TEMP1_AS_TOP_BC). Symétrique de predict_true_surface_signal, mais vers
    le bas : l'amplitude DIMINUE et la phase AUGMENTE avec la profondeur
    (l'inverse de l'extrapolation vers la surface) :

        A_deeper     = A_anchor * exp(-a_pred * dz_deeper)
        phase_deeper = phase_anchor + b_pred * dz_deeper
        T_deeper(t) = tendance_ref(t) + A_deeper*cos(w*t - phase_deeper) + résidu_ref(t)

    Motivation : si on réutilise le capteur le plus profond TEL QUEL comme
    condition limite du bas sans l'extrapoler, le point de comparaison
    immédiatement au-dessus (dans le schéma décalé, l'ancien capteur le plus
    profond sert à la fois de point de comparaison ET, dupliqué, de condition
    limite) devient rigoureusement identique à cette condition limite - le
    modèle ne peut alors le reproduire que si l'atténuation sur ce dernier
    intervalle est nulle, ce qui n'est jamais physiquement le cas et biaise le
    calage. La tendance lente et le résidu large-bande du capteur de référence
    sont conservés tels quels (une composante lente s'atténue très peu sur
    10cm supplémentaires, cf. discussion sur la condition limite haute).

    Parameters:
    - t_seconds: temps (s), commun à reference_values et à l'ajustement anchor_fit.
    - reference_values: série brute du capteur le plus profond (ex: Temp4).
    - anchor_fit: sortie de fit_sinusoid() pour ce même capteur.
    - dz_deeper: distance (m) à extrapoler au-delà de ce capteur.
    - a_pred, b_pred: taux d'atténuation/déphasage extrapolés à cette nouvelle
      profondeur (voir extrapolate_attenuation_phase).

    Returns:
    - dict avec 'predicted' (array, température prédite plus profond),
      'amplitude_pred', 'phase_pred', et 'amplitude_raw'/'phase_raw' (harmonique
      brute du capteur de référence, pour comparaison).
    """
    decomp = decompose_signal(t_seconds, reference_values, period_seconds, trend_order)
    w = decomp['w']
    t = np.asarray(t_seconds, dtype=float)

    amplitude_pred = anchor_fit['amplitude'] * np.exp(-a_pred * dz_deeper)
    phase_pred = anchor_fit['phase'] + b_pred * dz_deeper

    harmonic_pred = amplitude_pred * np.cos(w * t - phase_pred)
    predicted = decomp['trend'] + harmonic_pred + decomp['residual']

    return {
        'predicted': predicted,
        'amplitude_pred': amplitude_pred,
        'phase_pred': phase_pred,
        'amplitude_raw': decomp['amplitude'],
        'phase_raw': decomp['phase'],
    }


def depth_profile_regression(fits, depths, period_seconds, Cv, Cw=4.18e6):
    """
    Version "pédagogique" de l'analyse : au lieu de calculer les diffusivités
    apparentes paire par paire (ce qui multiplie les nombres et les combinaisons),
    on trace directement les deux droites théoriques prédites par le modèle de
    Stallman, en utilisant TOUTES les profondeurs à la fois :

        ln(amplitude(z)) = ln(A0) - a*z         -> une droite de pente -a
        phase(z)          = phase(0) + b*z       -> une droite de pente +b

    a et b sont donc simplement les pentes de deux régressions linéaires. C'est
    l'équivalent visuel direct de la théorie : plus la pente d'amplitude est
    forte, plus le signal s'atténue vite avec la profondeur (sol peu diffusif,
    ou écoulement contraire qui "pousse" la chaleur) ; la pente de phase donne le
    retard supplémentaire (en radians) pris par le signal à chaque mètre de
    profondeur supplémentaire.

    Parameters:
    - fits: dict {nom_profondeur: résultat de fit_sinusoid(...)}, incluant la
      profondeur 0 (référence de surface) si disponible.
    - depths: dict {nom_profondeur: profondeur en m, positif vers le bas}.
    - period_seconds: période utilisée pour les ajustements.
    - Cv, Cw: chaleurs volumiques du milieu saturé et de l'eau [J/m3/K].

    Returns:
    - dict contenant, entre autres, les tableaux z/log_amp/phase_unwrapped (pour
      tracer les points et les droites de régression), les pentes a et b, les
      coefficients de détermination R² (qualité de l'alignement des points sur
      la droite -> indicateur direct d'homogénéité du milieu), Gamma_ampl,
      Gamma_phase, et la vitesse de Darcy inversée uz.
    """
    names = sorted(depths, key=lambda n: depths[n])  # de la surface vers la profondeur
    z = np.array([depths[n] for n in names], dtype=float)
    amplitude = np.array([fits[n]['amplitude'] for n in names])
    phase_raw = np.array([fits[n]['phase'] for n in names])
    # Le signal profond est physiquement toujours en retard sur le signal de surface :
    # la phase doit croître de façon monotone avec z. np.unwrap corrige les sauts de 2*pi
    # dus au fait que arctan2 renvoie une valeur entre -pi et +pi.
    phase_unwrapped = np.unwrap(phase_raw)

    log_amp = np.log(amplitude)
    slope_amp, intercept_amp = np.polyfit(z, log_amp, 1)
    slope_phase, intercept_phase = np.polyfit(z, phase_unwrapped, 1)
    a = -slope_amp  # amortissement : ln(amp) DIMINUE avec z, d'où le signe -
    b = slope_phase  # déphasage : la phase AUGMENTE avec z

    def r_squared(x, y, slope, intercept):
        y_pred = slope * x + intercept
        ss_res = np.sum((y - y_pred) ** 2)
        ss_tot = np.sum((y - y.mean()) ** 2)
        return 1 - ss_res / ss_tot if ss_tot > 0 else np.nan

    r2_amp = r_squared(z, log_amp, slope_amp, intercept_amp)
    r2_phase = r_squared(z, phase_unwrapped, slope_phase, intercept_phase)

    w = 2 * np.pi / period_seconds
    gamma_ampl = w / (2 * a**2) if a > 0 else np.nan
    gamma_phase = w / (2 * b**2) if b > 0 else np.nan
    if a > 0 and b > 0:
        v, gamma_joint = invert_v_gamma(a, b, period_seconds)
        uz = darcy_velocity_from_v(v, Cv, Cw)
    else:
        v = gamma_joint = uz = np.nan
    uz_m_s = uz

    coherence = abs(gamma_ampl - gamma_phase) / gamma_phase if gamma_phase else np.nan

    return {
        'names': names, 'z': z, 'log_amp': log_amp, 'phase_unwrapped': phase_unwrapped,
        'slope_amp': slope_amp, 'intercept_amp': intercept_amp,
        'slope_phase': slope_phase, 'intercept_phase': intercept_phase,
        'a': a, 'b': b, 'r2_amp': r2_amp, 'r2_phase': r2_phase,
        'gamma_ampl_1e-6': gamma_ampl * 1e6, 'gamma_phase_1e-6': gamma_phase * 1e6,
        'coherence_rel_diff': coherence,
        'gamma_joint_1e-6': gamma_joint * 1e6 if np.isfinite(gamma_joint) else np.nan,
        'v_thermal_m_s': v,
        'uz_darcy_m_s': uz_m_s,
        'uz_darcy_mm_day': uz * 1000 * 86400 if np.isfinite(uz) else np.nan,
        'uz_darcy_mm_yr': uz * 1000 * 86400 * 365.25 if np.isfinite(uz) else np.nan,
    }
