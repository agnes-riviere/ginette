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




#######################
## INITIALIZE SCRIPT ##
#######################

# initialization: dos2unix in case on all files in folder
#find . -type f -exec dos2unix -q {} \;
cp E_zone_parameter_backup.dat E_zone_parameter.dat

## first thing check if values to replace are marked for replacement, otherwise throw an error
line1=$(sed -n '1p' < E_zone_parameter.dat) # zone 1
line2=$(sed -n '2p' < E_zone_parameter.dat) # zone 2
line3=$(sed -n '3p' < E_zone_parameter.dat) # zone 3
line4=$(sed -n '4p' < E_zone_parameter.dat) # zone 4
line5=$(sed -n '5p' < E_zone_parameter.dat) # zone 5



if [[ $line1 != *"["* ]] || [[ $line2 != *"["* ]] ||[[ $line3 != *"["* ]] ||[[ $line4 != *"["* ]] ||[[ $line5 != *"["* ]]
then
	echo "ERROR : SIMULATIONS CANNOT START :'( "
	echo "At least one parameter is not initialized to be replaced."
	exit 123
fi
gfortran -o traite_geom traite_geom.f90

# compile ginette
gfortran -o ginette ginette_V2.f


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

	# first check whether the simulation has already been run
	if [ ! -f ../S_piezoB_RD_${j}.dat ]; then

		# check whether simulations already converged
		if [ -e $fileNotConverging ] && grep -Fxq "$1_$j" $fileNotConverging ; then

			echo "jumping over line $j in file tested_values (didn't converge)"

		else # launch the simulations : not done and not too long

			# read parameters at line j
			echo "Reading line "$j
			paramSet=$(sed -n ${j}p < $file_i)
			echo "The parameters are " $paramSet

			# replace all  parameters in input files
			# cp E_zone_parameter.dat E_zone_parameter_backup.dat
			source replaceParams.sh $paramSet
			echo "..Physical parameters replaced"

			# delete Ginette output files in case
			# -f is force, no error if file doesnt exist
			rm -f S_temp*
			rm -f S_piezo*
			rm -f S_surf_piezo*
			# launch simulations for parameter set in line j
			# allow for max 15min
			 ./ginetteSteadyTransient_velocity.sh
				exitStatus=$?


				if [ $exitStatus -eq 124 ]; then # ginette had to stop before the end

					echo "Line $j didn't finish within time."
					echo "$1_$j" >> $fileNotConverging # store idx of file and idx of line in file
##.
				else # simulations finished, move files
					echo "Ginette exhausted"
#./OUTPUT/S_temp_
				NAMETP=""$OUT""$Rac"_"$TP"_"
#./OUTPUT/S_
				NAMEP=""$OUT""$Rac"_"
				echo $NAMEP
#./OBS/Obs_
				NAMEOT=""$OBS""$RacO"_"


				# rename output files

##.				mv S_temp_zns1.dat "$OUT"/S_temp_zns1_${j}.dat
##.				mv S_temp_zns2.dat "$OUT"/S_temp_zns1_${j}.dat
##.				mv S_temp_zns3.dat "$OUT"/S_temp_zns1_${j}.dat
##.				mv S_temp_TRes1.dat "$OUT"/S_temp_zns1_${j}.dat
##.				mv S_temp_TRes2.dat "$OUT"/S_temp_zns1_${j}.dat
##.				mv S_temp_HoboRD.dat "$OUT"/S_temp_zns1_${j}.dat
##.				mv S_temp_HoboRG.dat "$OUT"/S_temp_zns1_${j}.dat
					cp E_zone_parameter.dat "$PARAMETERS"/E_zone_parameter_${j}.dat
					mv S_piezoB_RD.dat "$OUT"/S_piezoB_RD_${j}.dat
					mv S_piezoB_RG.dat "$OUT"/S_piezoB_RG_${j}.dat
					mv S_surf_piezo.dat "$OUT"/S_surf_piezo_${j}.dat
					cp E_zone_parameter_backup.dat E_zone_parameter.dat
			for i in piezoB_RD piezoB_RG surf_piezo
			do
     		awk -f InGinette.awk -v  id_sim="$j" out="$NAMEP" obs="$NAMEOT" id="$i" anal_ginette.COMM
			mv awk.out anal_ginette.COMM
			./CmpKro/CmpKro0.01 anal_ginette.COMM apoub.log
			mv apoub.log ./SENSI/Sensi_"$j"_$i.dat
			done
			fi # close running simulation $1 $j

		fi # close check whether $1 $j had already gone over time

	else # the simulations have alread been run

		echo "line $j is already done!"

	fi

done

echo "file "$1" treated!"
