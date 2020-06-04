This User_guide will guide you step by step to help you complete the inversion based on the mini-LOMOS sensors, depending on how you chose to use it. 
Please refer to "LOMOS-mini: A coupled system quantifying transient water and heat exchanges in streambeds" from Cucchi et al. (2017) if you want to know more about it.
_____________________________________________________________________________________________________________________________

### Prelude :

Before making anything, please make sure you've installed the "dos2unix" UNIX package. You can do it from the terminal:
> apt install dos2unix

Also, make sure that all the files of the "mini-LOMOS" file are executable:
> chmod 775 *

In order to make the R scripts work, you'll also need the packages "lubridate", "hydroGOF" and "stringr".
You can respectively install them on R with the following command lines:
> install.packages(lubridate)

> install.packages(hydroGOF)

> install.packages(stringr)
_____________________________________________________________________________________________________________________________
### 1) Add field data.

There must be 2 ".csv" files (one for pressure and one for temperatures) with the following format: sensor-name_point_dd_MM_YY.csv
 
#### Particularly pay attention of the dates formatting in the titles, this could explain why the 1D model doesn't work.
-------------------------------------
example (respectively for temperature and pressure):
```
t502_Point1_07_11_14.csv
p502_Point1_07_11_14.csv
```
The format of the pressure differential file is:
#,dates, pressure_differential [m], temperature in stream [C]
```
1,04/11/2014 17:00:00,1.184,11.516,
etc.
```
The format of the temperature file is:
#,dates, temperature depth 1 [C], temperature depth 2 [C], temperature depth 3 [C], temperature depth 4 [C]
```
1,04/11/2014 17:00:00,12.775,12.92,13.69,13.112,
etc.
```
_____________________________________________________________________________________________________________________________
### 2) Complete the file : inversion.COMM
This file will be used by an R script to read the data.

Format the inversion.COMM command file columns as following :
```
 $1 : point name (no space) ex: Point1
 $2 : year (two last characters YY) ex for 2014: 14 
 $3 : month (two last characters MM) ex for November: 11 
 $4 : day (two last characters DD) ex: 07 
 $5 : name temperature sensor (four first characters) ex: t502
 $6 : name pressure sensor (four first characters) ex: p502
 $7 : time step (s) ex: 900
 $8 : nb of observations ex: 3 if the 4 PT100 work, 2 if one is broken, etc.
``` 
$8 = nb of working PT100-1; as the deepest working PT100 is used as a boundary condition. 
The upper boundary condition is the temperature sensor in the stream, in the pressure differential file.
```
 $9 : number of zones (clay, sand...) ex: 1
```
More than 2 would require more PT100 so it wouldn't make a lot of sense here.
```
 $10 : thickness of the upper zone (m)
```
$10 = distance between the riverbed and the bottom of the upper zones. Type in "0" if there is just 1 zone.
```
 $11: maximum duration of the inversion (s) ex : 864000 (= 10 days)
```
If the duration of the experiment is smaller than the prescribed duration, the duration of the experiment will be used (i.e. you can put 8640000000 to be sure to use the duration of the experiment).

------------------------------------
Example with every PT100 working and 2 zones:
```
Point1 14 11 07 t502 p504 900 3 2 0.20 864000
```
Example with 2 PT100 working and 1 zone:
```
Point43 17 05 18 t520 p520 900 2 1 0 864000
```
_____________________________________________________________________________________________________________________________
### 3) Complete the file : inversion_PT100.COMM

Each line contains the space between two following working PT100
The first line is the space between the riverbed and the first working PT100, the second between the first working PT100 and the second working PT100, etc.

------------------------------------
example for 4 PT100 with a space of 10 cm each:
```
10
10
10
10
```
example for 3 PT100 with a space of 10 cm each, but third PT100 is broken:
```
10
10
20
```
_____________________________________________________________________________________________________________________________
### 4) Complete the following file: inversion_parameter.COMM

Column 1 = parameters of interest:

k = intrinsic permeability [m2]

n = porosity 

l = solid thermal conductivity [W m−1 K−1]

r = solid density [kg m-3]

The bulk volumetric heat capacity of the porous medium is calculated  by the following equation :

**c_mr_m = c_w r_w n + c_s r (1-n)**

c_w = specific heat capacity of water [334 000 J kg−1 K−1]

r_w = water density [1 000 kg m-3] 

c_s = specific heat capacity of solid [J kg−1 K−1]

Column 2 = number of zones (e.g. clay, sand...)

Column 3 and 4 (respectively min and max) = range for the value test

Column 5 = number of tests within the previously defined range

#### Column 5 will be responsible for the number of simulations that will be launched.

-----------------------------------
example for 1 zone:
```
k 1 0001D-15 0001D-11 0003
n 1 0.40 0.40 0001
l 1 1300D-03 8400D-03 0002
r 1 2600D+00 2600D+00  0001
```
example for 2 zones:
```
k 1 0001D-15 0001D-11 0001
n 1 0.05 0.40 0001
l 1 1300D-03 8400D-03 0001
r 1 2600D+00 2600D+00  0001
k 2 0001D-15 0001D-11 0001
n 2 0.05 0.40 0001
l 2 1300D-03 8400D-03 0001
r 2 2600D+00 2600D+00  0001
```
_____________________________________________________________________________________________________________________________
### 5) Launch the screening of multiple experiments of temperature and pressure in the HZ with the awk script inversion.awk from the terminal:
> awk -f inversion.awk inversion.COMM
_____________________________________________________________________________________________________________________________
### 6) OUTPUT
The simulated results per mesh are in _ginette/application/mini-LOMOS/GINETTE_SENSI/OUTPUT/_

A comparative table of simulated versus observed is generated in _ginette/application/mini-LOMOS/SENSI/_ by "Comparaison_mailles_sim-obs.R". It is named "Results_Stats_sim-obs".

**Don't forget to copy/paste it in a specific repository if you want to keep track of these files, as launching another simulation will overwrite them.**
_____________________________________________________________________________________________________________________________
### 7) PLOTS
**Data check:**

The R script "formatToGinette.R" reads and interpolate all the data before they are used by ginette. 

At the end of the script (l.150), you'll find lines that will plot the interpolated data and save it in _ginette/application/mini-LOMOS/PLOT/Data_check_

You can do it before to use ginette to check if your data are correct, just remove the "#" using ctrl + shift + C and launch the script.

**Results:**

When you achieved satisfying calibrations, you can plot the simulated vs observed temperature time-series thanks to the temperature_time_series.R script in ginette/application/mini-LOMOS/PLOT/

It will also save the plots in _ginette/application/mini-LOMOS/PLOT/Data_check_
