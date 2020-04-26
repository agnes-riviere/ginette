# call by format_ts.sh

# $1 : sensor name (pressure or temperature)
# $2 : date d'enregistrement des donnees
# $3 : ouput file name ([test_nom].COMM du point ou f1.txt)

echo "-- Entering format_one_ts.sh"

nf=`awk -F"," -f clean.awk $1_$2.csv` #backquote = substitution de commande
echo $nf "fields in" $1"_"$2".csv" # execution de la commande
echo $nf > $3 # execution une 2e fois

for i in `seq 1 $nf`;
do
# pastes time series created in clean.awk to .dat files
mv out_$i.dat $1_$i.dat
echo $1_$i.dat > tmp.txt
cat $3 tmp.txt >tmp2.txt
mv tmp2.txt $3
done
rm tmp.txt

# file $3 contains the number of fields
# and the name of .dat files containing the series for $1_$2
