#!/bin/bash
# $1 : number of obervation
#program fortran qui genere la liste des param> tested_values
OUTPUT="OUTPUT"
TP="temperature"
#Lancement des simulation avec
chmod 755 *.sh
echo "Generating parameters' values"
gfortran -o values Calibrated_values.f95
./values
echo "Compiling gigi"
mkdir PARAMETERS
mkdir $OUTPUT
mkdir SENSI
echo "Configuring the links"
#rm CmpKro0.01
#ln -s ./CmpKro/CmpKro0.01


./main_oneFile.sh $2 $1
#./post_treat.sh tested_values $1
