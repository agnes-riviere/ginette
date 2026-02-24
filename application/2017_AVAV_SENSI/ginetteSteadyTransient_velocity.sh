# This file calls ginette in transient and steady state modes
OUT="./OUTPUT/"
OBS="./OBS/"
RacO="Obs"
Rac="Sim"
TP="temperature"
P="charge"

# compile ginette
gfortran -o AvAv ginette_V2.f

# change parameters in input files for steady state
# echo "Ginette's steady"

# option -i remplace directement le fichier donne en entree
# option -e permet de passer plusieurs commandes a la suite
# sed -i -e "s/rp=1/rp=0/g" E_parametre.dat
# sed -i -e "s/rpth=1/rpth=0/g" E_p_therm.dat
# sed -i -e "s/ichi2=1/ichi2=0/g" E_cdt_initiale.dat
# sed -i -e "s/itempi=1/itempi=0/g" E_cdt_initiale.dat
# sed -i -e "s/iclchgt=1/iclchgt=0/g" E_cdt_aux_limites.dat

# call ginette steady state
##./ginette_velocity

echo "Ginette transient"

# change parameters in input files for transient
#sed -i -e "s/rp=0/rp=1/g" E_parametre.dat
##sed -i -e "s/rpth=0/rpth=1/g"  E_p_therm.dat
#sed -i -e "s/ichi=0/ichi=1/g" E_cdt_initiale.dat
##sed -i -e "s/itempi=0/itempi=1/g" E_cdt_initiale.dat
#sed -i -e "s/iclchgt=0/iclchgt=1/g" E_cdt_aux_limites.dat

# paste pressure results from steady state
# awk '{print $4}' S_pression_charge_temperature.dat > E_pression_initiale.dat

# call ginette transient
#./traite_geom

./AvAv
echo "Ginette exhausted"



