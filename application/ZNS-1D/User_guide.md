## GINETTE User guide
AUCUNE MODIFICATION DU CODE N'IMPLIQUE LES DEVELLOPEURS AGNES RIVIERE JULIO GONCALVES ANNE JOST

_Authors contacts_
- agnes.riviere@mines_paristech.fr
- goncalves@cerege.fr
- anne.jost@upmc.fr

Pour toutes utilisation de ce code, les références suivantes doivent être citées:

> _Rivière, A., Jost, A., Gonçalvès, J., and Font, M. (2019) – Experimental and numerical studies of pore-water pressure variations in sub-permafrost groundwater. Cold Regions Science and Technology, 158, 76-94. ⟨10.1016/j.coldregions.2018.11.005⟩. _

> _Rivière, A., Gonçalvès, J., Jost, A. and Font, M. (2014) – Experimental and numerical assessment of stream-aquifer exchanges during disconnection. Journal of Hydrology, 517, 574-583. ⟨10.1016/j.jhydrol.2014.05.040⟩._

### FICHIER ENTREE
    A REMPLIR :
    E_cdt_aux_limites.dat
    E_cdt_initiale.dat
    E_parametre.dat
    E_p_therm.dat
    E_geom.dat

    OPTIONAIRE
    E_charge_initiale.dat : charge initiale par colonne
    E_pression_initiale.dat : 1 pression intiale par maille
    E_temperature_initiale.dat : temperature initiale par maille

#### To compile 
    gfortran ../../src/ginette_V2.f -o ginette
#### to run
    ./ginette
#### PARAMETER AND SET UP  
#### Flow model 
#### E_parametre.dat
    dt=000900D+00	s	time sted ; format d9.0
    nitt=0000000001		final time = nitt*unite
    unite=8640D+01    s    unite
    l=1000D+00	m 	lenght    ; format d8.0
    az=00011D+00   depth      ;   format d9.0
    imaille=0	      =0 Read mesh: E_voisins.dat, E_coordonnees.dat, E_def_maille.dat  imaille=1 automatic building
    akx=0256D-13	m2	horizontal intrinsic permeability
    akz=0256D-13	m2	vertical intrinsic permeability
    omp=0.3800	     total porosity

    recorded time step itsortie*unitsortie
    itsortie=00086400
    unitsortie=00001D+00
    ytest=ZNS      	cas test TH2 TH3 THL ANU DTS AVA ZHR ZHZ ZNS1D
#### INITIAL CONDITIONS  
#### E_cdt_initiale.dat
    chgi=00117D+00 	m		homogeneous hydraulic head
    tempi=0012d+00	C		homogeneous temperature
    ichi=1			ichi=1 hydrostatic hydraulic head E_charge_initiale.dat one colunm with the value (Pa) of hydraulic head
    ichi2=0				ichi2=1 spatial variation of the pressure E_pression_initiale.dat one value by celle
    itempi=0		itempi=1 spatial variation of the temperature E_temperature_initiale.dat, one temp ny cell
#### BOUNDARY CONDITIONS  
##### E_cdt_aux_limites
    Homogenous :
    Left:
    icl_gauche=-1		icl=-1 neuman (fluxes), icl=-2 dirichlet (hydraulic head)
    valcl_gauche=00000000d+00	units : fluxes (m/s) , hydraulic head (m)
    Right
    icl_droite=-2
    valcl_droite=00000065d-02
    Top
    icl_haut=-1
    valcl_haut=00000000d+00
    icl_bas=-1
    Bottom
    valcl_bas=00000000d+00
    iclchgt=1		special conditions in the subroutine variation_cdt_limites in ginette_V2.f

    TEST ZNS
    fichier input
    E_p_therm.dat : fichier thermique et nom du test
    E_parametre.dat : parametre et set up



    Conditions initiales
    E_cdt_initiale.dat
    E_charge_initiale.dat


    Conditions limites
    E_cdt_aux_limites.dat


### OUTPUT  
    S_hydraulic_conductivities_profil_t.dat'     time, z, K (m/s) 
	S_saturation_profil_t.dat                    time, z , water saturation
    S_pressure_profil_t.dat						 time, z, pressure (Pa)










