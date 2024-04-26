This User_guide will guide you step by step to help you complete the inversion based on the mini-LOMOS sensors, depending on how you chose to use it. 
Please refer to "LOMOS-mini: A coupled system quantifying transient water and heat exchanges in streambeds" from Cucchi et al. (2017) if you want to know more about it.
_____________________________________________________________________________________________________________________________

1) Complete the following file: inversion_parameter.COMM

Column = parameters of interest:
k= intrinsic permeability
n= porosity
l= heat capacity
r= medium density

Column 2 = number of zones (e.g. clay, sand...)
Column 3 and 4 (respectively min and max) = range for the value test (litterature)
Column 5 = number of tests within the previously defined range

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

2) Complete the file : inversion.COMM

Format the inversion.COMM command file columns as following :

 $1 : point name (no space) ex: Point1
 $2 : year (two last characters YY) ex for 2014: 14 
 $3 : month (two last characters MM) ex for November: 11 
 $4 : day (two last characters DD) ex: 07 
 $5 : name temperature sensor (four first characters) ex: t502
 $6 : name pressure sensor (four first characters) ex: p502
 $7 : thickness of the studied hyporheic zone (= distance between the 2 furthest working PT100, in meters) ex: 0.40
 $8 : time step (s) ex: 900 
 $9 : nb of meshes between the river bed and the second PT100 ex: 11 (= 11 cm with a mesh of 1 cm)
 $10 : nb of meshes between the river bed and the third PT100 ex:21
 $11 : nb of meshes between the river bed and the fourth PT100 ex:31
 $12 : nb of functioning PT100 (3 if every sensor works, 2 if one is broken, etc.)
 $13 : number of zones (clay, sand...). More than 2 would require more PT100 so it wouldn't make a lot of sense here)
 $14 : thickness of the upper zone (= distance between the riverbed and the bottom of the upper zone, in meters)

------------------------------------
example with every PT100 working and 2 zones :
Point1 14 11 07 t502 p504 0.40 900 11 21 31 3 2 0.20
_____________________________________________________________________________________________________________________________

3) Complete the file : inversion_PT100.COMM

Each line contains the space between two following PT100
The first line is the space between the riverbed and PT100 1, the second between PT100 1 and PT100 2, etc.

------------------------------------
example for 4 PT100 with a space of 10 cm each:
10
10
10
10
_____________________________________________________________________________________________________________________________

4) Add field data with the following format: sensor-name_point_dd_MM_YY.csv 

-------------------------------------
example :
t502_Point1_07_11_14.csv

The format of the pressure differential file is :
#,dates,pressure_differential [m],temperature in stream [C],
1,04/11/2014 17:00:00,1.184,11.516,

The format of the temperature file is :
#,dates,temperature depth 1 [C],temperature depth 2 [C],temperature depth 3 [C],temperature depth 4 [C],
1,04/11/2014 17:00:00,12.775,12.92,13.69,13.112,
_____________________________________________________________________________________________________________________________

5) Launch the screening of multiple experiments of temperature and pressure in the HZ with the awk script inversion.awk from the terminal:
> awk -f inversion.awk inversion.COMM

