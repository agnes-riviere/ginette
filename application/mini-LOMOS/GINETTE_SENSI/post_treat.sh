#!/bin/bash

for i in agglo_out.awk remove_first_line.sh remove_first_line.awk agglo_line.sh cat_agglo.sh $1
do
cp $i ./SENSI
done
cd SENSI
rm tmp*
rm awk_agglo.out
echo "Now agglomerating sensi estimates in " `pwd`
echo "Check everything's here " `ls *.awk *.sh`
awk -f agglo_out.awk -v nsensors=$2 $1
sed -i -e 's/,/ /g' awk_agglo.out
cat ../header.in awk_agglo.out > Sensi.dat
cat ../header.values tmp.out > Value.dat
join -1 1 -2 1 Value.dat Sensi.dat > ../Sensi_final.dat
cd ..
sed -i -e 's/ /,/g' Sensi_final.dat
echo "Sensibility values extracted in Sensi_final.dat"

#Il faudra tout effacer quand le script marchera completement. A Valider maintenant
#rm -rf ./SENSI
#rm -rf ./OUTPUT
rm -rf ./PARAMETERS
#rm -rf ./OBS
