#!/bin/bash
# this file runs Ginette with parameter sets given in $1 th line of parameterSets.txt
# it handles one set of porosity - heat conductivity - heat capacity
# it loops over values of permeability provided in parameters/permeabilitySets/

OUT="./OUTPUT/"
OBS="./OBS/"
RacO="Obs"
Rac="Sim"
TP="temperature"
P="charge"
nbzone=$1
var1=1


#######################
## INITIALIZE SCRIPT ##
#######################

# initialization: dos2unix in case on all files in folder
#find . -type f -exec dos2unix -q {} \;
			if [ $nbzone -eq $var1 ]
    		then
			cp E_zone_parameter_backup.dat E_zone_parameter.dat
			else
			cp E_zone_parameter_backup_2zones.dat E_zone_parameter.dat
			fi



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
	if [ ! -f ../Sim_temperature_j${j}_maille1_t.dat ]; then

		# check whether simulations already converged
		if [ -e $fileNotConverging ] && grep -Fxq "$1_$j" $fileNotConverging ; then

			echo "jumping over line $j in file tested_values (didn't converge)"

		else # launch the simulations : not done and not too long

			# read parameters at line j
			echo "Reading line "$j
			paramSet=$(sed -n ${j}p < $file_i)
			echo "The parameters are " $paramSet  $1

			if [ $nbzone -eq $var1 ]
    		then
			cp E_zone_parameter_backup.dat E_zone_parameter.dat
			else
			cp E_zone_parameter_backup_2zones.dat E_zone_parameter.dat
			fi

			if [ $nbzone -eq $var1 ]
    				then
			source replaceParams.sh $paramSet
			echo "..Physical parameters replaced"
			else
			source replaceParams_2zones.sh $paramSet
			echo "..Physical parameters replaced for 2 zones"
			fi

			# delete Ginette output files in case
			# -f is force, no error if file doesnt exist
			#rm -f Sim_temperature_maille*

			# launch simulations for parameter set in line j
			# allow for max 15min
			./ginetteSteadyTransient_velocity.sh
			exitStatus=$?


			if [ $exitStatus -eq 124 ]; then # ginette had to stop before the end

				echo "Line $j didn't finish within time."
				echo "$1_$j" >> $fileNotConverging # store idx of file and idx of line in file

			else # simulations finished, move files
				echo "Ginette exhausted"

				NAMETP=""$OUT""$Rac"_"$TP"_maille"
#./OUTPUT/Sim_temperature_maille
				NAMEP=""$OUT""$Rac"_"$P"_maille"
#./OUTPUT/Sim_charge_maille
#Supprime les espaces multiples et les espaces en debut et fin de chaque ligne
			sed "{s/\ \ */\ /g;s/^\ *//g;s/ $//g}" Sim_temperature_maille1_t.dat > Sim_temperature_maille1_t_clean.dat
			sed "{s/\ \ */\ /g;s/^\ *//g;s/ $//g}" Sim_temperature_maille2_t.dat > Sim_temperature_maille2_t_clean.dat
			sed "{s/\ \ */\ /g;s/^\ *//g;s/ $//g}" Sim_temperature_maille3_t.dat > Sim_temperature_maille3_t_clean.dat
			rm Sim_temperature_maille1_t.dat
			rm Sim_temperature_maille2_t.dat
			rm Sim_temperature_maille3_t.dat
			mv Sim_temperature_maille1_t_clean.dat "$OUT"/Sim_temperature_maille1_${j}.dat
	        	mv Sim_temperature_maille2_t_clean.dat "$OUT"/Sim_temperature_maille2_${j}.dat
		        mv Sim_temperature_maille3_t_clean.dat "$OUT"/Sim_temperature_maille3_${j}.dat
		        mv S_vitesse_nmaille2_hb.dat "$OUT"/S_vitesse_nmaille2_hb_${j}.dat
	        	mv Sim_velocity_profil_t.dat "$OUT"/Sim_velocity_profil_t_${j}.dat
	        	mv Sim_heat_flux_profil_t.dat "$OUT"/Sim_heat_flux_profil_t_${j}.dat
	        	mv Sim_temperature_profil_t.dat "$OUT"/Sim_temperature_profil_t_${j}.dat



#			nb_obs=$2
# 			for i in `seq 1 $nb_obs`
# 			do
# 			awk -f InGinette.awk -v  id_sim="$j" out="$NAMETP" obs="$NAMEOT" id="$i" anal_ginette.COMM
# 			mv awk.out anal_ginette.COMM
# 			./CmpKro/CmpKro0.01 anal_ginette.COMM apoub.log
# 			mv apoub.log ./SENSI/Sensi_temperature_"$j"_$i.dat
# 			done
			fi # close running simulation $1 $j


		fi # close check whether $1 $j had already gone over time

	else # the simulations have alread been run

		echo "line $j is already done!"

	fi

done

echo "file "$1" treated!"
