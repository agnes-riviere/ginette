# -*- coding: utf-8 -*-
"""
Solutions analytiques DIRECTES (forward) utilisées pour valider Ginette contre
la théorie, indépendamment de tout site réel — voir application/heat_transport_analytical_validation/.

Contrairement à Stallman_analysis.py (qui INVERSE des profils de température
mesurés pour estimer v/Gamma), ce module part de paramètres physiques CONNUS
(v, Gamma, ou log_k/lam/poro) et calcule la température T(z,t) ou T(z) prédite
par la théorie, pour comparaison directe à la sortie numérique de Ginette.

Trois solutions, correspondant aux trois cas de application/heat_transport_analytical_validation/ :

1. Stallman (1965) — régime TRANSITOIRE, signal de surface sinusoïdal,
   milieu homogène. Relation directe (inverse mathématique de
   Stallman_analysis.invert_v_gamma) :

       a + i*b = (-v + sqrt(v**2 + 4*i*Gamma*w)) / (2*Gamma)
       T(z,t) = T0 + A0*exp(-a*z)*cos(w*t - b*z)

2. Bredehoeft & Papadopulos (1965) — régime PERMANENT, milieu homogène,
   températures fixées aux deux bords (Dirichlet haut ET bas) :

       T(z) = T0 + (TB-T0) * (exp(q*z/gamma) - 1) / (exp(q*L/gamma) - 1)

3. Shan & Bodvarsson (2004) — extension MULTI-COUCHES de la précédente
   (Bredehoeft & Papadopulos en est le cas particulier à 1 couche) :

       Ti(z) = Ci1*exp(q*z/gamma_i) + Ci2   (z dans la couche i)

   Implémentée ici dans la version RE-PRÉSENTÉE ET CORRIGÉE par Kurylyk,
   Irvine, Carey, Briggs, Werkema & Bonham (2017) "Heat as a groundwater
   tracer in shallow and deep heterogeneous media: Analytical solution,
   spreadsheet tool, and field applications", Hydrological Processes 31,
   2648-2661 (eq. 7-12 du papier - mêmes coefficients Ci1/Ci2/alpha que
   Shan & Bodvarsson 2004, mais gamma corrigé, voir plus bas, et nombre de
   Péclet thermique ajouté, voir thermal_peclet_number()). Kurylyk et al.
   valident cette solution par comparaison à SUTRA (éléments finis) sur des
   systèmes à 1/2/3/4 couches - même principe de validation croisée que
   application/heat_transport_analytical_validation/ pour Ginette.

Convention commune aux cas 2 et 3 (reprise directement des deux papiers) :
- q : flux de Darcy [m/s], POSITIF vers le bas (recharge). z est la profondeur,
  POSITIVE vers le bas, z=0 au sommet du domaine (PAS la convention "altitude"
  z_top=0/z_bottom<0 utilisée ailleurs dans ce projet - convertir avant appel :
  z_profondeur = z_top - z_altitude).
- gamma_i = lambda_i / Cw (PAS la diffusivité thermique réelle lambda_i/Cv :
  Shan & Bodvarsson (2004) définissent gamma avec la capacité calorifique de
  L'EAU seule, Cw, pas celle du milieu saturé Cv - point de confusion signalé
  explicitement par Kurylyk et al. (2017), corrigé ici en conséquence).
- Kurylyk et al. (2017), eq. 12 : le nombre de Péclet thermique d'un système à
  n couches (rapport flux advectif/flux conductif) utilise comme conductivité
  "bulk" de l'ensemble la moyenne HARMONIQUE des conductivités de couche
  pondérée par leur épaisseur ("because the heat flow is perpendicular to the
  layering", citant Lunardini 1981) - EXACTEMENT le même principe physique
  que le fix de matt() dans src/ginette_V2.f90 (2026-07-23, voir
  project_matt_conductivity_interface_fix dans la mémoire du projet) :
  confirmation indépendante, par la littérature spécifique au traçage
  thermique en milieu stratifié, que la moyenne harmonique inter-couches est
  la bonne pratique physique à l'interface entre deux conductivités
  différentes.
"""
import numpy as np
from scipy.special import erf, erfc
from scipy.optimize import least_squares

# Constantes reprises EXACTEMENT des valeurs codées en dur dans le template
# Ginette E_p_therm_bck.dat (alandae, cpe, cpm - voir
# GINETTE_SENSI/E_p_therm_bck.dat de 1D_Stream_aquifer_GridSearch), pour que
# les prédictions analytiques utilisent la même physique que ce que Ginette
# calcule réellement en interne, plutôt que des valeurs génériques de la
# littérature qui pourraient légèrement différer.
CW_VOL = 4.185e6     # capacité calorifique volumique de l'eau [J/m3/K] (cpe=4185 * rho_eau=1000)
KW = 0.598           # conductivité thermique de l'eau [W/m/K] (alandae=0598D-03 kg.m/s-3/C)
CS_SPECIFIC = 1000.0  # capacité calorifique massique du solide [J/kg/K] (cpm=00001000, GLOBAL,
                       # pas calibrée par zone - seuls k1/n1/l1/r1 varient dans le grid search)


def bulk_thermal_conductivity(ks, poro, kw=KW):
    """
    Conductivité thermique du milieu SATURÉ (k_bulk), à partir de la
    conductivité de la fraction SOLIDE seule (ks, le paramètre `lam` calibré
    dans le grid search) et de la porosité - loi de mélange ARITHMÉTIQUE,
    celle utilisée par défaut par Ginette (ymoycondtherm=ARITH, voir
    Stallman_analysis.solid_thermal_conductivity dont c'est l'inverse exact) :

        k_bulk = kw*poro + ks*(1-poro)
    """
    return kw * poro + ks * (1 - poro)


def volumetric_heat_capacity(poro, ref_density, cs_specific=CS_SPECIFIC, Cw_vol=CW_VOL):
    """
    Capacité calorifique volumique du milieu SATURÉ (Cv, [J/m3/K]), moyenne
    arithmétique pondérée par la porosité entre l'eau et la fraction solide
    (ref_density = REF_DENSITY dans les scripts Ginette, densité du grain
    minéral en kg/m3, ex: 2650 pour du quartz) :

        Cv = poro*Cw_vol + (1-poro)*ref_density*cs_specific
    """
    return poro * Cw_vol + (1 - poro) * ref_density * cs_specific


def darcy_velocity_from_head(log_k, dh, L, temperature_celsius=15.0, rho=1000.0, g=9.81):
    """
    Vitesse de Darcy uz [m/s] produite par une différence de charge dh [m]
    imposée sur une longueur L [m], pour une perméabilité intrinsèque
    k=10**log_k [m2] - relation directe (inverse de
    Stallman_analysis.intrinsic_permeability) :

        K = k*rho*g/mu(T)   ;   uz = K*dh/L

    Sert à choisir les conditions limites hydrauliques (h_top, h_bottom) de
    Ginette pour obtenir un flux de Darcy CONNU q=uz, nécessaire pour prédire
    T(z) avec bredehoeft_papadopulos_profile / shan_bodvarsson_layered_profile.
    """
    mu = 2.414e-5 * 10 ** (247.8 / (temperature_celsius + 133.15))  # viscosité de l'eau (Vogel)
    k = 10 ** log_k
    K = k * rho * g / mu
    return K * dh / L


def stallman_forward(v, gamma, period_seconds):
    """
    Relation DIRECTE de Stallman (1965) : amortissement a (1/m) et déphasage
    b (rad/m) prédits à partir de la vitesse frontale thermique v (m/s,
    = uz*Cw/Cv) et de la diffusivité thermique gamma (m2/s, = k_bulk/Cv) -
    inverse mathématique exact de Stallman_analysis.invert_v_gamma :

        a + i*b = (-v + sqrt(v**2 + 4*i*gamma*w)) / (2*gamma)

    Returns: (a, b), tous deux positifs (v=0 -> cas purement diffusif, a=b).
    """
    w = 2 * np.pi / period_seconds
    lam = (-v + np.sqrt(v**2 + 4j * gamma * w)) / (2 * gamma)
    return lam.real, lam.imag


def shan_bodvarsson_layered_profile(z, q, layer_thicknesses, layer_conductivities, T0, TB, Cw_vol=CW_VOL):
    """
    Solution de Shan & Bodvarsson (2004), régime PERMANENT, milieu à n
    couches de conductivités thermiques distinctes (Bredehoeft & Papadopulos
    1965 en est le cas particulier n=1 - voir docstring du module et
    bredehoeft_papadopulos_profile ci-dessous, qui appelle directement cette
    fonction avec une seule couche).

    Parameters:
    - z: profondeur(s) (m, POSITIVE vers le bas, 0 au sommet du domaine) où
      évaluer T - scalaire ou array.
    - q: flux de Darcy (m/s), POSITIF vers le bas (recharge).
    - layer_thicknesses: épaisseurs des couches (m), de la plus proche de la
      surface à la plus profonde (liste/array de longueur n).
    - layer_conductivities: conductivité thermique BULK (k_bulk, W/m/K -
      utiliser bulk_thermal_conductivity() si on part de lam/poro calibrés
      Ginette) de chaque couche, même ordre que layer_thicknesses.
    - T0: température au sommet du domaine (z=0) [°C].
    - TB: température au bas du domaine (z=somme des épaisseurs) [°C].

    Returns:
    - T(z), même forme que l'array z en entrée.
    """
    b = np.asarray(layer_thicknesses, dtype=float)
    lam = np.asarray(layer_conductivities, dtype=float)
    n = len(b)
    gamma = lam / Cw_vol  # eq. 2 de Shan & Bodvarsson (2004) - PAS la diffusivité réelle

    d = np.cumsum(b)  # profondeur cumulée jusqu'au bas de chaque couche (d[-1] = L)
    d_top = np.concatenate(([0.0], d[:-1]))  # profondeur du haut de chaque couche

    if q == 0:
        # Cas purement conductif (q->0) : eq. 7-11 dégénèrent (division par q dans
        # alpha), la solution limite est le profil linéaire par morceaux classique
        # (flux conductif constant, continuité de température ET de flux aux
        # interfaces - cas particulier bien connu, pas besoin de alpha/C1/C2).
        flux_conductif = (TB - T0) / np.sum(b / lam)
        T_interfaces = np.concatenate(([T0], T0 + np.cumsum(b / lam) * flux_conductif))
        z_arr = np.atleast_1d(np.asarray(z, dtype=float))
        T = np.empty_like(z_arr)
        for i in range(n):
            mask = (z_arr >= d_top[i]) & (z_arr <= d[i] + 1e-9)
            T[mask] = T_interfaces[i] + (z_arr[mask] - d_top[i]) / lam[i] * flux_conductif
        return T if hasattr(z, "__len__") else T[0]

    alpha = np.exp(q * np.sum(b / gamma))  # eq. 11

    C2 = (alpha * T0 - TB) / (alpha - 1)  # eq. 8 (identique pour toutes les couches)
    C1 = np.zeros(n)
    C1[0] = (TB - T0) / (alpha - 1)  # eq. 9
    for i in range(n - 1):
        C1[i + 1] = np.exp(q * d[i] * (1.0 / gamma[i] - 1.0 / gamma[i + 1])) * C1[i]  # eq. 10

    z_arr = np.atleast_1d(np.asarray(z, dtype=float))
    T = np.empty_like(z_arr)
    for i in range(n):
        mask = (z_arr >= d_top[i]) & (z_arr <= d[i] + 1e-9)
        T[mask] = C1[i] * np.exp(q * z_arr[mask] / gamma[i]) + C2
    return T if hasattr(z, "__len__") else T[0]


def thermal_peclet_number(q, layer_thicknesses, layer_conductivities, Cw_vol=CW_VOL):
    """
    Nombre de Péclet thermique d'un système à n couches (Kurylyk et al. 2017,
    eq. 12) - rapport du flux de chaleur advectif au flux conductif :

        Pe = q * sum(bi / gamma_i)     avec gamma_i = lambda_i / Cw_vol

    Équivalent à ln(alpha) (voir eq. 11 dans shan_bodvarsson_layered_profile).
    Se réduit à q*L/gamma pour n=1 (cas Bredehoeft & Papadopulos 1965).
    Note : Kurylyk et al. expriment aussi Pe de façon équivalente via la
    conductivité "bulk" du système, calculée comme la moyenne HARMONIQUE des
    conductivités de couche pondérée par leur épaisseur (voir docstring du
    module) - non utilisée ici, la formule directe en sum(bi/gamma_i) est
    plus simple et strictement équivalente.
    """
    b = np.asarray(layer_thicknesses, dtype=float)
    lam = np.asarray(layer_conductivities, dtype=float)
    gamma = lam / Cw_vol
    return q * np.sum(b / gamma)


def lunardini_freezing_profile(x, t, Tm, T0=4.0, Ts=-6.0, Tf=0.0,
                                k1=3.4644, k2=2.9414, k3=2.4184,
                                C=690360.0, gammad=1680.0, xif=0.0782, xi0=0.2000,
                                Lf=334720.0):
    """
    Solution analytique de Lunardini (1988, CRREL Report 82-2) pour la
    propagation d'un front de gel en régime TRANSITOIRE dans un milieu
    poreux semi-infini initialement non gelé - cas dit "exact" (valeurs
    non physiques choisies pour coller exactement aux hypothèses
    mathématiques de Lunardini - conductivités k1/k2/k3 et capacité
    volumique C CONSTANTE dans les 3 zones, voir docstring de
    ymoycondtherm="LUNAR"/ytest="THL" dans ginette_V2.f90, qui implémentent
    directement ce cas particulier).

    3 zones (x = profondeur depuis la surface, positive vers le bas) :
    - zone 1 (gelée)      : 0 <= x <= X1,   T1(x) = Ts + (Tm-Ts)*erf(psi*x/X1)/erf(psi)
    - zone 2 ('mushy')    : X1 <= x <= X2,  cf. ci-dessous (erf ou erfc selon gamma)
    - zone 3 (non gelée)  : x >= X2,        T3(x) = T0 - (T0-Tf)*erfc(beta*x/X2)/erfc(beta)

    ATTENTION au signe de la 1ère équation du système non-linéaire : c'est
    bien un MOINS devant a1*c1*erf(psi)/(...). Certaines versions publiées
    de cette solution (dont McKenzie et al. 2007) y mettent un PLUS, ce qui
    ne redonne pas les valeurs de référence. Vérifié en reproduisant
    psi/gamma/beta/X1/X2 pour Tm=-4°C et Tm=-1°C à 1e-5 près avec le signe
    MOINS ci-dessous - la 2nde équation, elle, est correcte telle quelle.

    Parameters:
    - x: profondeur(s) [m], positive vers le bas, où évaluer T - scalaire ou array.
    - t: temps [s] depuis le changement brutal de température en surface.
    - Tm: température solidus (limite bas de la zone mushy - haut de la zone
      gelée) [°C] - PARAMÈTRE VARIABLE du cas test (Tableau 1/2/4 : -0.1,
      -1.0 ou -4.0°C). Correspond à tsg/tsd dans E_p_therm.dat.
    - T0, Ts, Tf: températures initiale (loin de la surface), imposée en
      surface, et liquidus (haut de la zone mushy = 0°C par convention
      Lunardini) [°C]. Tf correspond à tlg/tld dans E_p_therm.dat.
    - k1, k2, k3: conductivités thermiques BULK des 3 zones [W/m/K]
      (gelée/mushy/non gelée) - valeurs par défaut = celles hardcodées dans
      ginette_V2.f90 pour ymoycondtherm="LUNAR".
    - C: capacité calorifique volumique (identique dans les 3 zones, avant
      ajout du terme de chaleur latente apparente dans la zone mushy)
      [J/m3/K] - valeur par défaut = celle hardcodée dans ginette_V2.f90
      pour ytest="THL" (690360).
    - gammad: densité sèche du solide [kg/m3].
    - xif, xi0: ratio massique eau non gelée/solide sec, zones gelée et non
      gelée respectivement (xif < xi0).
    - Lf: chaleur latente de fusion PAR KG D'EAU [J/kg] (convention
      Lunardini - PAS la même que le paramètre "lat" de E_p_therm.dat, qui
      attend la chaleur latente par kg de GLACE : lat = Lf*rho_eau/rho_glace,
      soit 334720 J/kg (eau) -> 363826 J/kg (glace) pour
      rho_glace=920 kg/m3).

    Returns:
    - (T, X1, X2): T de même forme que x, et les profondeurs du front gel/mushy
      (X1) et mushy/non-gelé (X2) [m] à l'instant t.
    """
    alpha1 = k1 / C
    alpha3 = k3 / C
    alpha4 = k2 / (C + gammad * Lf * (xif - xi0) / (Tm - Tf))

    a1 = np.sqrt(alpha1 / alpha4)
    a2 = np.sqrt(alpha4 / alpha3)
    b1 = (Tm - Ts) / (Tf - Tm)
    b2 = (Tm - Tf) / (T0 - Tf)
    c1 = k2 / k1
    c2 = k3 / k2

    def eqs(v):
        psi, gamma = v
        F1 = b1 * np.exp(-(psi**2) * (1 - a1**2)) - a1 * c1 * erf(psi) / (erf(gamma) - erf(a1 * psi))
        F2 = b2 * np.exp(-(gamma**2) * (1 - a2**2)) + a2 * c2 * (erf(gamma) - erf(a1 * psi)) / (1 - erf(a2 * gamma))
        return [F1, F2]

    best = None
    for psi0 in (0.02, 0.05, 0.1, 0.2, 0.3, 0.5):
        for gamma0 in (0.3, 0.6, 1.0, 1.5, 2.0, 3.0):
            res = least_squares(eqs, [psi0, gamma0], bounds=([1e-8, 1e-8], [15, 15]), xtol=1e-14, ftol=1e-14)
            if best is None or res.cost < best.cost:
                best = res
    psi, gamma = best.x
    beta = gamma * np.sqrt(alpha4 / alpha3)

    X1 = psi * np.sqrt(4 * alpha1 * t)
    X2 = gamma * np.sqrt(4 * alpha4 * t)

    x_arr = np.atleast_1d(np.asarray(x, dtype=float))
    T = np.empty_like(x_arr)

    frozen = x_arr <= X1
    thawed = x_arr >= X2
    mushy = ~frozen & ~thawed

    T[frozen] = Ts + (Tm - Ts) * erf(psi * x_arr[frozen] / X1) / erf(psi)
    T[thawed] = T0 - (T0 - Tf) * erfc(beta * x_arr[thawed] / X2) / erfc(beta)
    if gamma > 1:
        T[mushy] = Tf + (Tm - Tf) * (erfc(gamma * x_arr[mushy] / X2) - erfc(gamma)) / \
            (erfc(gamma * X1 / X2) - erfc(gamma))
    else:
        T[mushy] = Tf + (Tm - Tf) * (erf(gamma * x_arr[mushy] / X2) - erf(gamma)) / \
            (erf(gamma * X1 / X2) - erf(gamma))

    T = T if hasattr(x, "__len__") else T[0]
    return T, X1, X2


def bredehoeft_papadopulos_profile(z, q, k_bulk, L, T0, TB, Cw_vol=CW_VOL):
    """
    Solution de Bredehoeft & Papadopulos (1965), régime PERMANENT, milieu
    HOMOGÈNE - cas particulier à une seule couche de
    shan_bodvarsson_layered_profile (voir docstring du module) :

        T(z) = T0 + (TB-T0) * (exp(q*z/gamma) - 1) / (exp(q*L/gamma) - 1)
        gamma = k_bulk / Cw_vol

    Parameters:
    - z: profondeur(s) (m, positive vers le bas, 0 au sommet) où évaluer T.
    - q: flux de Darcy (m/s), positif vers le bas (recharge).
    - k_bulk: conductivité thermique du milieu SATURÉ (W/m/K) - utiliser
      bulk_thermal_conductivity() si on part de lam/poro calibrés Ginette.
    - L: épaisseur totale du domaine (m).
    - T0, TB: températures fixées en haut (z=0) et en bas (z=L) [°C].
    """
    return shan_bodvarsson_layered_profile(z, q, [L], [k_bulk], T0, TB, Cw_vol=Cw_vol)
