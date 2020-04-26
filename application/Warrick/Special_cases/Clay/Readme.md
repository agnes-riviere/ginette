## GINETTE User guide
AUCUNE MODIFICATION DU CODE N'IMPLIQUE LES DEVELLOPEURS AGNES RIVIERE JULIO GONCALVES ANNE JOST

_Authors contacts_for ginette code
- agnes.riviere@mines_paristech.fr
- goncalves@cerege.fr
- anne.jost@upmc.fr

_Contributor to this test:
- deniz.kilic@mines-paristech.fr

Pour toutes utilisation de ce code, les références suivantes doivent être citées:

> _Rivière, A., Jost, A., Gonçalvès, J., and Font, M. (2019) – Experimental and numerical studies of pore-water pressure variations in sub-permafrost groundwater. Cold Regions Science and Technology, 158, 76-94. ⟨10.1016/j.coldregions.2018.11.005⟩. _

> _Rivière, A., Gonçalvès, J., Jost, A. and Font, M. (2014) – Experimental and numerical assessment of stream-aquifer exchanges during disconnection. Journal of Hydrology, 517, 574-583. ⟨10.1016/j.jhydrol.2014.05.040⟩._



## TEST WAR
Warrick, A. W., Lomen, D. O., & Yates, S. R. (1985). A generalized solution to infiltration. Soil Science Society of America Journal, 49(1), 34-38.
    fichier input

#### E_parametre.dat
dt=30s
Ksref (m/s)	=1.00E-06
n_VG= 1.53 --> dz max=0.20m
alpha_vg (1/m)	= 1
Residual saturation = 0.5139


Column thickness= 1 meter
dz=0025D-04	m


#### E_cdt_aux_limites.dat
Prescribed recharge
icl_haut=-1													c Surface	C ECOULEMENT
valcl_haut=00000452d-07	valeur flux (m/s) valeur charge (m)								c Surface  	C ECOULEMENT
