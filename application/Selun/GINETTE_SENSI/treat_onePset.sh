#!/bin/bash
OUT="./OUTPUT/"
OBS="./OBS/"
RacO="Obs"
Rac="Sim"
TP="temperature"
P="charge"


awk -f traite_ginette.awk -v n=$1 Ss=$2 kh=$3 kv=$4 E_parametre.dat
awk -f traite_ginette_therm.awk -v lambdas=$5 cs=$6 rhos=$7 E_p_therm.dat
cp awk.ginette ./PARAMETERS/E_parametre\_$1\_$2\_$3\_$4\_$5\_$6\_$7.dat
mv  awk.ginette E_parametre.dat
cp awk.ginette_therm ./PARAMETERS/E_p_therm\_$1\_$2\_$3\_$4\_$5\_$6\_$7.dat
mv  awk.ginette_therm E_p_therm.dat
echo "Ginette's steady"
sed -i -e "s/rp=1/rp=0/g" E_parametre.dat
sed -i -e "s/rpth=1/rpth=0/g" E_p_therm.dat
sed -i -e "s/ichi2=1/ichi2=0/g" E_cdt_initiale.dat
sed -i -e "s/itempi=1/itempi=0/g" E_cdt_initiale.dat
sed -i -e "s/iclchgt=1/iclchgt=0/g" E_cdt_aux_limites.dat
./ginette
echo "steady state done"
echo "Ginette's steady"
awk  '{print $4}' S_pression_charge_temperature.dat>E_pression_initiale.dat
#echo "Ginette's steady"awk  '{print $6}' S_pression_charge_temperature.dat > E_temperature_initiale.dat
sed -i -e "s/rp=0/rp=1/g" E_parametre.dat
sed -i -e "s/rpth=0/rpth=1/g"  E_p_therm.dat
sed -i -e "s/ichi2=0/ichi2=1/g" E_cdt_initiale.dat
sed -i -e "s/itempi=0/itempi=1/g" E_cdt_initiale.dat
sed -i -e "s/iclchgt=0/iclchgt=1/g" E_cdt_aux_limites.dat
echo "Ginette transient"
./ginette
echo "Ginette exhausted"
NAMETP=""$OUT""$Rac"_"$TP"_maille"
NAMEP=""$OUT""$Rac"_"$P"_maille"
NAMEOT=""$OBS""$RacO"_"$TP"_maille"

#cp S_pression_charge_temperature.dat $OUT\/S\_$1\_$2\_$3\_$4\_$5\_$6\_$7.dat


echo "Analysing Ginette's performances with RAC" $NAMEP $NAMEOT
for i in 1 2 3
do
#cp Sim_charge_maille$i\_t.dat "$NAMEP"$i\_$1\_$2\_$3\_$4\_$5\_$6\_$7.dat
echo "cp Sim_temperature_maille"$i"_t.dat "$NAMETP""$i"_"$1"_"$2"_"$3"_"$4"_"$5"_"$6"_"$7".dat"
cp Sim_temperature_maille$i\_t.dat "$NAMETP"$i\_$1\_$2\_$3\_$4\_$5\_$6\_$7.dat
dos2unix "$NAMETP"$i\_$1\_$2\_$3\_$4\_$5\_$6\_$7.dat
awk -f InGinette.awk -v a="$1" b="$2" c="$3" d="$4" e="$5" f="$6" g="$7" out="$NAMETP" obs="$NAMEOT" id="$i" anal_ginette.COMM
mv awk.out anal_ginette.COMM
./CmpKro/CmpKro0.01 anal_ginette.COMM apoub.log
mv apoub.log ./SENSI/$TP$i\_$1\_$2\_$3\_$4\_$5\_$6\_$7.dat
#rm "$NAMETP"$i\_$1\_$2\_$3\_$4\_$5\_$6\_$7.dat
done
