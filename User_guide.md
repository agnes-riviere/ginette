## GINETTE User guide
AUCUNE MODIFICATION DU CODE N'IMPLIQUE LES DEVELLOPEURS AGNES RIVIERE JULIO GONCALVES ANNE JOST

_Authors contacts_
- agnes.riviere@minesparis.psl.eu
- nicolas.radic@minesparis.psl.eu
- goncalves@cerege.fr
- anne.jost@upmc.fr

Pour toutes utilisation de ce code, les références suivantes doivent être citées:

> _Rivière, A., Jost, A., Gonçalvès, J., and Font, M. (2019) – Experimental and numerical studies of pore-water pressure variations in sub-permafrost groundwater. Cold Regions Science and Technology, 158, 76-94. ⟨10.1016/j.coldregions.2018.11.005⟩. _

> _Rivière, A., Gonçalvès, J., Jost, A. and Font, M. (2014) – Experimental and numerical assessment of stream-aquifer exchanges during disconnection. Journal of Hydrology, 517, 574-583. ⟨10.1016/j.jhydrol.2014.05.040⟩._

### COMPILATION
#### Use gfortan
$ gorftran -o #NAME_EXECUTABLE# ginette_V2.f90

#### use meson 
sudo apt-get install python3 ninja-build meson
$ cd /path/to/ginette/src
$ meson setup builddir && cd builddir
$ meson compile
$ meson test


### FICHIER ENTREE
    A REMPLIR :
    E_cdt_aux_limites.dat
    E_cdt_initiale.dat
    E_parametre.dat
    E_p_therm.dat
    E_geom.dat

#### OPTIONAIRE
    E_charge_initiale.dat : charge initiale par colonne
    E_pression_initiale.dat : 1 pression intiale par maille
    E_temperature_initiale.dat : temperature initiale par maille

#### maillage manuel :
    E_row.dat le numero de ligne de chaque maille
    E_colonne.dat le numero de ligne de chaque maille
    E_voisins.dat
    E_coordonnee.dat
    E_def_maille.dat
    E_BordRD.dat            identifiant mailles bord droite
    E_BordRG.dat             identifiant mailles bord gauche
    zonage de parametres
    E_zone.dat  numero de zone de chaque maille
    E_zone_parameter.dat     parametre de chaque zone

#### Conditions limites heterogene
    E_cl_ecoulement.dat type BC Droite gauche haut bas valeurs BC Droite gauche haut bas
    E_cl_thermique.dat

    E_BordRD.dat            identifiant mailles bord droite
    E_BordRG.dat             identifiant mailles bord gauche
    E_Id_river.dat          identifiant maille riviere
     E_Id_river_max.dat     identifiant maille maximum riviere

    E_chargeT_RD.dat       valeur charge des mailles E_BordRD en fonction du temps
    E_chargeT_RG.dat       valeur charge des mailles E_BordRG en fonction du temps
    E_chargeT_Riv.dat      valeur charge des mailles E_BordRD en fonction du temps

#### To compile 
    gfortran ginette_V2.f -o ginette
#### to run
    ./ginette
#### create model 
    gfortran traite_geom.f90 -o treat_geom
    ./treat_geom
#### to debug
    gfortran -g -fbacktrace -ffpe-trap=zero,overflow,underflow ginette_V2.f -o ginette
    gdb ./ginette
    run
#### PARAMETER AND SET UP  
#### Flow model 
#### E_parametre.dat
    ec=1                   iec=1 flow active
    rp=1			rp=1 transient state ; rp=0 steady state
    dt=000900D+00	s	time sted ; format d9.0
    nitt=0000000001		final time = nitt*unite
    unite=8640D+01    s    unite
    l=1000D+00	m 	lenght    ; format d8.0
    az=00011D+00   depth      ;   format d9.0
    imaille=0	      =0 Read mesh: E_voisins.dat, E_coordonnees.dat, E_def_maille.dat  imaille=1 automatic building
    itopo=1          itopo=1 mesh in function of E_geom.dat
    akx=0256D-13	m2	horizontal intrinsic permeability
    akz=0256D-13	m2	vertical intrinsic permeability
    omp=0.3800	     total porosity
    ss=3.80D-01	     storage coefficient/porosity water model unconfined without unsatured process
    itsortie=00086400    recorded time step itsortie*unitsorti
    unitsortie=00001D+00
    solver=CGS         BIC BiConjugate Gradient, CGS Conjugate Gradient Squared Method, LIB inria library (n'est pas active actuellement, pb lorsque ecoulement + transport chaleur)


### unsaturated flow
#### E_parametre.dat
omp=0.3630		porosity                                                             
yunconfined=UNS 
ivg=1	
ans=0153D-02		n Van Genuchten parameter
asp=0100D-02	1/m	alpha Van Genuchten paramerer 
swres=0.5124		residual saturation
Be careful
- dz should respect this criteria dz<= 1/(5*alpha) in the capilarity fring dans la frange capillaire
- dt=dz/vfront with vfront the saturation front velocoty vfront=infiltration/(porosité*(smax-smin)

### Heat model 
#### E_parametre.dat
    ith=1				ith=0 heat inactive ith=1 heat active
    time sted
    dt=000900D+00	s	pas de temps ; format d9.0
    final time = nitt*unite
    nitt=0000000001		nombre d'iterations en temps
    unite=8640D+01    s    unite iteration temps
    l=1000D+00	m 	lenght    ; format d8.0
    az=00011D+00   depth      ;   format d9.0
    imaille=0	      =0 Read mesh: E_voisins.dat, E_coordonnees.dat, E_def_maille.dat  imaille=1 automatic building
    itopo=1          itopo=1 mesh in function of E_geom.dat
    recorded time step itsortie*unitsortie
    itsortie=00086400
    unitsortie=00001D+00
    solver=CGS         BIC BiConjugate Gradient, CGS Conjugate Gradient Squared Method, LIB inria library
#### E_p_therm.dat
    ithec=1               ithec=1 advection +conduction ithec=0 conduction without advection warning iec in E_parametre.dat should be equal to 0
    alandam=7.70D+00 	  solid thermal conductivity
    rhosi=2650D+00	kg/m3  solid density
    cpm=0835D+00	m2/s2/C	specific heat capacity
    ytest=RIV      	cas test TH2 TH3 THL ANU DTS AVA ZHR ZHZ
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
