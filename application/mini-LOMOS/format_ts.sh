# appele par format_ginette.sh
# $1 : point name
# $2 : temperature sensor name
# $3 : beginning of the experiment, the date must be formatted as dd_mm_yy
# $4 : pressure sensor name
# $5 :  beginning of the experiment, the date must be formatted as dd_mm_yy

echo "-- Entering format_ts.sh"
./format_one_ts.sh $4_$1 $5 test_$1.COMM
./format_one_ts.sh $2_$1 $3 f1.txt

# fusion des fichiers f1.txt et test_[nomDuPoint].COMM
cat f1.txt test_$1.COMM > f2.txt
mv f2.txt test_$1.COMM
rm f1.txt

# reste fichier test_[nom du point].COMM
# il contient
# - nb de fields en temperature
# - une ligne correspondant au nom du fichier .dat * nb de fields
# - nb de fields en pression
# - une ligne correspondant au nom du fichier .dat * nb de fields
