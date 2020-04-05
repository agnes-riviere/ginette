The scripts work with libts version 1.68  and with libpc version 1.24. Watch the variables V_TS V_LP in the Makefiles

1) You MUST to complete the file : inversion_parameter.COMM
parameters zone min max nb of tested values
n=  : min :max  nb of tested value
k=   : min max intrinsic permeability nb of tested values
l= :min thermal conductivity nb of tested values
r=  : min max heat capicity  nb of tested values
-----------------------------------
example for 2 zones :
k 1 0001D-15 0001D-11 0001
n 1 0.05 0.40 0001
l 1 1300D-03 8400D-03 0001
r 1 2600D+00 2600D+00  0001
k 2 0001D-15 0001D-11 0001
n 2 0.05 0.40 0001
l 2 1300D-03 8400D-03 0001
r 2 2600D+00 2600D+00  0001

_____________________________________________________________________________________________________________________________
2) Complete inversion.COMM
You MUST format your command file as inversion.COMM with this following format
# $1 : name point ex:Point1
# $2 : year (two last characters YY) ex: 14
# $3 : month (two last characters MM) ex:  11
# $4 : day (two last characters DD) ex: 07
# $5 : name temperature file (four first characters) ex: t502
# $6 : name pressure sensor (four first characters) ex: p504
# $7 : thickness of the studied hyporheic zone (m)ex: 0.40
# $8 : duration of one day in second ex:86400
# $9 : delta t (s) ex: 900
# $10 : nmaille1 PT100 1 ex: 11
# $11 : nmaille2 PT100 2 ex:21
# $12 : nmaille3 PT100 3 ex:31
# $13 : Number of observation PT100
# $14 : Number of area 1 or 2
# $15 : Limit between the two area

------------------------------------
example:
Point1 14 11 07 t502 p504 0.40 86400 900 11 21 31 3 2 0.20

3) Complete the file
inversion_PT100.COMM
Each line contains the space between two PT100
The first line is the space between the river and the PT100 1 and the last one between the deepest recorded PT100 and the other one.
------------------------------------
example for 4 PT100 with a space of 10 cm:
10
10
10
10
_____________________________________________________________________________________________________________________________

4) Add the field data of the point the name must be :
name sensor_name point_day_month_year.csv example t502_Point1_07_11_14.csv
The format of the pressure differential file is :
-------------------------------------
#,dates,pressure_differential [m],temperature in stream [C],
1,04/11/2014 17:00:00,1.184,11.516,
The format of the temperature file is :
#,dates,temperature depth 1 [C],temperature depth 2 [C],temperature depth 3 [C],temperature depth 4 [C],
1,04/11/2014 17:00:00,12.775,12.92,13.69,13.112,


3) Launch the screening of multiple experiment of temperature and pressure in the HZ with the awk script inversion.awk:
> awk -f inversion.awk inversion.COMM

