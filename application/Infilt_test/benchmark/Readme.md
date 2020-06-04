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



## TEST
Explain the goal of the test!!!

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


### Benchmark test

# E_cdt_aux_limites.dat
icl_haut=-1	 c Surface
valcl_haut=00000000d-00	valeur flux (m/s) valeur charge (m)
icl_bas=-2	
valcl_bas=-0000008d+00	valeur flux (m/s) valeur charge (m)
itlecture=00086400
# E_cdt_initale.dat
chgi=-0008D+00 	m	charge initiale ; d9.0 
conci=0015d+00		concentration initiale ; d8.0			
tempi=0005D+00  C	temperature initiale ; d8.0
ichi=1			ichi=1 --> E_charge_initiale.dat
ichi2=0			ichi2=1 --> E_pression_initiale.dat

# E_charge_initiale.dat
- Steady state profile result from:
> initial condition = -8m
> top boundary condition = 0m
> bottom boundary condition = -8m

# E_charge_bas_t.dat
-8 m

# E_debit_haut_t.dat
3.96372e-09 m/s

# E_parametre.dat

dt=00900D+00	s
nitt=0000000001		nombre d'iterations en temps               
unite=3162D+04    s    unite iteration temps             al=0001D+00	m 	longueur
az=00070D+00   m  
ixz=1
imaille=1	      =0 lecture maillage : voisins, coordonnees, def_maille : largeur (am), hauteur (bm), nb colonnes (nci), nb lignes (nri);  =1 construction maillage automatique                                             
itopo=0                 itopo=1 variation de altitude de top bot par colonne pas de E_geom.dat
ibuilt=1			construction maillage manuelle et zonage des parametres
reptop=00000D+00 m	altitude top
repbot=-0070D+00 m	altitude bottom

akz=0589D-14	m2	permeabilite intrinseque verticale      
omp=0.2930		porosite
ans=0110D-02		n parametre de Van Genuchten
asp=0100D-02	1/m	alpha parametre de Van Genuchten
swres=0.0000		saturation residuelle



