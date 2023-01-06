#!/bin/bash
# All the following informations come from the file "inversion.COMM" and are imported here by the file "inversion.awk":

# $1 : point name
# $2 : beginning of the experiment
# $3 : temperature sensor name
# $4 : pressure sensor name
# $5 : time step (s)
# $6 : nb of functioning PT100 (3 if every sensor works, 2 if one is broken, etc.)
# $7 : Number of areas (clay, sand...)
# $8 : Thickness of each area (separator = "", in meters)

echo "-- Entering format_ginette.sh"
## Read calibrated parameter
input_parameter="inversion_parameter.COMM"
#fichier cree qu'on remplit


tested_ranges="E_tested_ranges.dat"
nbligne=$(wc -l $input_parameter |cut -d" " -f1)
a=$((${7}*4));

if [[ -f "$tested_ranges" ]]; then
    rm awk.ginette
rm E_tested_ranges.dat
fi
if [[ -f E_nb_zone.dat ]]; then
rm E_nb_zone.dat
fi

for i in `seq 1 $a`; do echo 1 1 1 >> awk.ginette; done

echo ${7}  ${8} >> E_nb_zone.dat

if  [ $nbligne != $a ]
then
        echo "Le nombre de zone et de parametres des fichiers de inversion_parameter COMM ne correspondent pas !"
        echo  " ${8} zones  sont imposees" "inversion_parameter COMM contient" $nbligne "au lieu de" $a
exit
else
        echo "Le nombre de zone et de parametres des fichiers de inversion_parameter COMM correspondent !" $a
fi


cp awk.ginette $tested_ranges
while read  param zone min max pas
do
nbzone=${7}
k="k"
n="n"
l="l"
r="r"

if [ $param = $k ]
then
awk -f traite_k.awk -v zone=$zone min=$min max=$max nbpas=$pas nbzone=$nbzone $tested_ranges
cp awk.ginette $tested_ranges
fi
if [ $param = $n ]
then
awk -f traite_n.awk -v zone=$zone min=$min max=$max nbpas=$pas nbzone=$nbzone $tested_ranges
cp awk.ginette $tested_ranges
fi
if [ $param = $l ]
then
awk -f traite_l.awk -v zone=$zone min=$min max=$max nbpas=$pas nbzone=$nbzone $tested_ranges
cp awk.ginette $tested_ranges
fi
if [ $param = $r ]
then
awk -f traite_r.awk -v zone=$zone min=$min max=$max nbpas=$pas nbzone=$nbzone $tested_ranges
cp awk.ginette $tested_ranges
fi
done <$input_parameter

chmod 755 *.sh
#supprime tous les fichiers contenant "test_" dans le repertoire ./ (non recursif), ils feraient planter "formatToGinette.R"
#find -maxdepth 1 -name '*test_*.COMM' | xargs rm
OBS=./GINETTE_SENSI/OBS/
if [ -d $OBS ]; then
rm -r $OBS
fi
mkdir $OBS
OUTPUT=./GINETTE_SENSI/OUTPUT
if [ -d $OUTPUT ]; then
rm -r $OUTPUT
fi
mkdir $OUTPUT

cp $tested_ranges ./GINETTE_SENSI/




## fichier [test_nomdupoint.COMM] cree
## Creation des fichiers d'entree pour Ginette
R CMD BATCH formatToGinette.R

##fichier observation pour inversion

##for i in `seq 1 ${6}`;
##do
##mv Obs_temperature_maille$i\_t.dat ./GINETTE_SENSI/OBS/Obs_temperature_maille$i\_t.dat
##done #done i

for i in E_pression_initiale.dat E_charge_t.dat E_temp_t.dat E_temperature_initiale.dat nitt.dat model.dat
do
mv $i ./GINETTE_SENSI/
done

mv E_nb_zone.dat ./GINETTE_SENSI/

cd ./GINETTE_SENSI
echo "Now preparing ginette in" `pwd`
dos2unix --quiet *.dat



nitt=`awk '{print $NF}' nitt.dat`
awk -f CmdConfig.awk -v nitt="$nitt" E_parametre.dat
mv awk.out E_parametre.dat



#Add the number of cell and size of the colunm
nm=`awk '{print $1}' model.dat`
L=`awk '{print $2}' model.dat`
nmaille1=`awk '{print $3}' model.dat`

if (( ${6} > 1 )) ; then
nmaille2=`awk '{print $4}' model.dat`
else
nmaille2=10
fi
if (( "${6}" > 2 )) ; then
nmaille3=`awk '{print $5}' model.dat`
else
nmaille3=10
fi

awk -f CmdConfig_nm.awk -v nm="$nm" -v  L="$L" E_parametre.dat
mv awk.out E_parametre.dat
awk -f CmdConfig_nmaille.awk -v nmaille1="$nmaille1" -v  nmaille2="$nmaille2"  -v  nmaille3="$nmaille3" E_parametre.dat
mv awk.out E_parametre.dat

echo "Launching the sensitivity for" $1
./HZ1D.sh ${6} ${7}
cd ..
# Calcul des critères stats
if [ -d SENSI ]; then
rm -r SENSI
fi
mkdir SENSI
R CMD BATCH Comparaison_mailles_sim-obs.R

echo "file "$1" treated!"
