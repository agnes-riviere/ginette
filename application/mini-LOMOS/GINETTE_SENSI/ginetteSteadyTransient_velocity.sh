#!/bin/bash
# This file calls ginette in transient and steady state modes
OUT="./OUTPUT/"
OBS="./OBS/"
RacO="Obs"
Rac="Sim"
TP="temperature"
P="charge"

# compile ginette
gfortran -o ginette_velocity ../../../src/ginette_V2.f 


echo pwd
pwd


echo "Ginette transient"

# change parameters in input files for transient
sed -i -e "s/rp=0/rp=1/g" E_parametre.dat
sed -i -e "s/rpth=0/rpth=1/g"  E_p_therm.dat
sed -i -e "s/ichi2=0/ichi2=1/g" E_cdt_initiale.dat
sed -i -e "s/itempi=0/itempi=1/g" E_cdt_initiale.dat
sed -i -e "s/iclchgt=0/iclchgt=1/g" E_cdt_aux_limites.dat

# paste pressure results from steady state
awk '{print $4}' S_pression_charge_temperature.dat > E_pression_initiale.dat

# call ginette transient
./ginette_velocity

echo "Ginette exhausted"



