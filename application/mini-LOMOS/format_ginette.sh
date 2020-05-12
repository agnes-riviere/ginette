#!/bin/bash
# All the following informations come from the file "inversion.COMM" and are imported here by the file "inversion.awk":

# $1 : point name
# $2 : temperature sensor name
# $3 : beginning of the experiment
# $4 : pressure sensor name
# $5 : beginning of the experiment
# $6 : time step (s)
# $7 : thickness of the studied hyporheic zone (= distance between the 2 furthest working PT100, in meters)
# $8 : nb of meshes between the river bed and the second PT100
# $9 : nb of meshes between the river bed and the third PT100
# $10 : nb of meshes between the river bed and the fourth PT100
# $11 : nb of functioning PT100 (3 if every sensor works, 2 if one is broken, etc.)
# $12 : Number of areas (clay, sand...)
# $13 : Thickness of each area (separator = "", in meters)

echo "-- Entering format_ginette.sh"
## Read calibrated parameter
input_parameter="inversion_parameter.COMM"
#fichier cree qu'on remplit
tested_ranges="E_tested_ranges.dat"
nbligne=$(wc -l $input_parameter |cut -d" " -f1)
a=$((${12}*4));
rm awk.ginette
rm E_tested_ranges.dat
for i in `seq 1 $a`; do echo 1 1 1 >> awk.ginette; done
rm E_nb_zone.dat
echo ${12}  ${13} >> E_nb_zone.dat

if  [ $nbligne != $a ]
then
        echo "Le nombre de zone et de parametres des fichiers de inversion_parameter COMM ne correspondent pas !"
        echo  " ${12} zones  sont imposees" "inversion_parameter COMM contient" $nbligne "au lieu de" $a
exit
else
        echo "Le nombre de zone et de parametres des fichiers de inversion_parameter COMM correspondent !" $a
fi


cp awk.ginette $tested_ranges
while read  param zone min max pas
do
nbzone=${12}
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
find -maxdepth 1 -name '*test_*.COMM' | xargs rm
mkdir ./GINETTE_SENSI/OBS/
rm -rf ./GINETTE_SENSI/OBS/*
rm -rf ./GINETTE_SENSI/OUTPUT/*
cp $tested_ranges ./GINETTE_SENSI/

## creation du fichier [test_nomdupoint.COMM]
./format_ts.sh $1 $2 $3 $4 $5
#time observation and simulation -> Comparision to adapt the boundary condition files
echo $6 > tmp.txt
cat tmp.txt test_$1.COMM > tmp2.txt
mv tmp2.txt test_$1.COMM
rm tmp.txt



#Add the number of cell for the interpolation of the initial value in the file main.c
L=$7
nm=$(awk -v L=$L 'BEGIN { print ((L ) / 0.01) }');
echo $nm > tmp.txt
cat tmp.txt test_$1.COMM > tmp2.txt
mv tmp2.txt test_$1.COMM
rm tmp.txt


echo $5 > tmp.txt
cat tmp.txt test_$1.COMM > tmp2.txt
mv tmp2.txt test_$1.COMM
rm tmp.txt

## fichier [test_nomdupoint.COMM] cree
## Creation des fichiers d'entree pour Ginette
R CMD BATCH formatToGinette.R

##fichier observation pour inversion
for i in `seq 1 ${11}`;
do
mv Obs_temperature_maille$i\_t.dat ./GINETTE_SENSI/OBS/Obs_temperature_maille$i\_t.dat
done #done i

for i in E_pression_initiale.dat E_charge_t.dat E_temp_t.dat E_temperature_initiale.dat nitt.dat
do
mv $i ./GINETTE_SENSI/
done

mv E_nb_zone.dat ./GINETTE_SENSI/

cd ./GINETTE_SENSI
echo "Now preparing ginette in" `pwd`
dos2unix *.dat
#echo "compiling CmpKro"
#cd ./CmpKro
#make clean
#make
#cd ../


nitt=`awk '{print $NF}' nitt.dat`
awk -f CmdConfig.awk -v nitt="$nitt" E_parametre.dat
mv awk.out E_parametre.dat



pwd
echo "Lenght of 1D model L=" $8
awk -f CmdConfig_nm.awk -v nm="$nm" -v  L="$L" E_parametre.dat
mv awk.out E_parametre.dat
awk -f CmdConfig_nmaille.awk -v nmaille1="${8}" -v  nmaille2="${9}"  -v  nmaille3="${10}" E_parametre.dat
mv awk.out E_parametre.dat

#echo "Launching the sensitivity for" $1
#./HZ1D.sh ${11} ${12}
cd ..
# Calcul des critères stats
if [ -d SENSI ]; then
rm -r SENSI
fi
#mkdir SENSI
#R CMD BATCH Comparaison_mailles_sim-obs.R

# rm -rf Sim_*.dat
#cp ./Sensi_final.dat ../$1_Sensi_final.dat
#echo "Sensibility file moved in "$1"_Sensi_final.dat"
