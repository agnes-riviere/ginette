# Launch the script with gawk -F"," -f clean.awk p504_Point1_07_11_14.csv

# Appele par format_one_ts.sh
# un seul argument : le nom du fichier de donnees csv

# creates file out_i.dat containing time series of data per field in csv
# prints of this file go in text file (third argument of format_one_ts.sh)

# output : number of fields in csv file

function isfloat(x)
{
return x ~ /^-?(([0-9]+\.[0-9]*)|(\.?[0-9]+))([eE][-+][0-9]+)?$/
}

BEGIN{  
    name = "out";
}
{
    if (isfloat($3)){
	for(i=3;i<NF;i++){
	    printf("%s %f\n",$2,$i)>name"_"(i-2)".dat";
		echo name"_"(i-2)".dat";
	    nfields=NF-3;
	}
    }

}
END{
    printf("%d\n",nfields);
}
