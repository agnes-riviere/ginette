## GINETTE User guide
AUCUNE MODIFICATION DU CODE N'IMPLIQUE LES DEVELLOPEURS AGNES RIVIERE JULIO GONCALVES ANNE JOST

_Authors contacts_
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
    E_p_therm.dat : fichier thermique et nom du test WAR
    E_parametre.dat : parametre et set up
#### To plot solution 
	R script

### INPUT FILE
    E_cdt_aux_limites.dat
    E_parametre.dat
    E_charge_initiale.dat


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
    az=00011D+00   depth      ;   format d9.0 positive value
	reptop=00000D-04 m	altitude top
	repbot=-00012+00 m	altitude bottom
	dz=0001D-02	m	discretisation
    akx=0256D-13	m2	horizontal intrinsic permeability
    akz=0256D-13	m2	vertical intrinsic permeability
    omp=0.3800	     total porosity
	ans=0203D-02		n parametre de Van Genuchten						CCCCCCCCCCCCCCCCCCCCCCCCCCC ZNS
	asp=0664D-02	1/m	alpha parametre de Van Genuchten 					CCCCCCCCCCCCCCCCCCCCCCCCCCC ZNS
	swres=0.0005		saturation residuelle							CCCCCCCCCCCCCCCCCCCCCCCCCCC ZNS
    recorded time step itsortie*unitsortie
    itsortie=00086400
    unitsortie=00001D+00
    ytest=ZNS      	cas test TH2 TH3 THL ANU DTS AVA ZHR ZHZ ZNS1D
#### INITIAL CONDITIONS  
#### E_charge_initiale.dat
one presure head by cell (P/rho/g)

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

### OUTPUT  
    S_hydraulic_conductivities_profil_t.dat'     time, z, K (m/s) 
	S_saturation_profil_t.dat                    time, z , water saturation
    S_pressure_profil_t.dat						 time, z, pressure (Pa), hydraulic head (m)










