# this file runs Ginette with parameter sets given in $1 th line of parameterSets.txt
# it handles one set of porosity - heat conductivity - heat capacity
# it loops over values of permeability provided in parameters/permeabilitySets/

OUT="./OUTPUT/"
OBS="./OBS/"
RacO="Obs"
Rac="S"
TP="temp"
Chg="charge"
Piez="piezo"
PARAMETERS="./PARAMETERS/"





#####################################
## INITIALIZE PARAMETERS AND FILES ##
#####################################

# store name of file where parameters are stored
file_i=tested_values

echo "..Considering "$file_i

# get number of parameter sets to try (ie number of lines in $file_i)
nbLines=$(cat $file_i | wc -l)

echo "File "$file_i" contains "$nbLines" lines (ie parameter sets)."

# initialize name of file where non converging parameter sets will be stored
fileNotConverging="../nonConvergingSets.txt"

#################
## ACTUAL LOOP ##
#################

for j in $(eval echo "{1..$nbLines}"); # loop over lines in $file_i
do


			# read parameters at line j
			echo "Reading line "$j
			paramSet=$(sed -n ${j}p < $file_i)
			echo "The parameters are " $paramSet

#./OUTPUT/S_temp_
				NAMETP=""$OUT""$Rac"_"$TP"_"
#./OUTPUT/S_
				NAMEP=""$OUT""$Rac"_"
				echo $NAMEP
#./OBS/Obs_
				NAMEOT=""$OBS""$RacO"_"
			for i in piezoB_RD piezoB_RG surf_piezo
			do
     		awk -f InGinette.awk -v  id_sim="$j" out="$NAMEP" obs="$NAMEOT" id="$i" anal_ginette.COMM
			mv awk.out anal_ginette.COMM
			./CmpKro/CmpKro0.01 anal_ginette.COMM apoub.log
			mv apoub.log ./SENSI/Sensi_"$j"_$i.dat
			done


done

echo "file "$1" treated!"
