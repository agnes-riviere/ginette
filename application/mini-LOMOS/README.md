Ginette application
=====================================================

# Authors:
*  Riviere , Agnes, agnes.riviere@mines_paristech.fr
*  Cucchi , karina, karina.cucchi@gmail.com 

# Contributions :
*  Schuller , Guillaume, gschuler@engees.eu
*  Flipo , Nicolas, nicolas.flipo@mines_paristech.fr


# created in 2008_. Copyright ©
- Sorbonne Universite
- Universite Aix Marseille
- Mines ParisTech

# References:
* Cucchi, K., Rivière, A., Flipo, N., A. Baudin, A. Berrhouma, F. Rejiba & Rubin, Y. (2018). LOMOS-mini: a coupled pressure and temperature system for local estimations of water and heat exchanges and sediment properties in streambeds. Journal of hydrology, 561, 1037-104. ⟨10.1016/j.jhydrol.2017.10.074⟩.
* Rivière, A., Gonçalvès, J., Jost, A., & Font, M. (2014). Experimental and numerical assessment of transient stream-aquifer exchange during disconnection. Journal of Hydrology, 517, 574–583. https://doi.org/10.1016/j.jhydrol.2014.05.040
* Rivière, A., Jost, A., Gonçalvès, J. & Font, M. (2018). Pore water pressure evolution below a freezing front under saturated conditions: Large-scale laboratory experiment and numerical investigation. Cold Regions Science and Technology, 158, 76-94. https://doi.org/10.1016/j.coldregions.2018.11.005
* Rivière, A., Gonçalvès, J., Jost, A., (2020). agnes-riviere/ginette: Ginette-2020-09 (Version 2020-09). Zenodo. http://doi.org/10.5281/zenodo.4058821

_**Update coming up later**_

# Inverting Temperature Hyporheic Zone Data for Water and Heat Flux Estimation

Using temperature data from the hyporheic zone in conjunction with pressure differential data, our scripts can estimate the water and heat flux that occurs between a surface stream and the subsurface. Python users and batch processing users can benefit from the script in their own unique ways. One user guide is for Python users and the other is for batch users. In the batch user guide, users are instructed on how to run the script via a command line interface, whereas in the Python user guide, users are instructed on how to use the script. Depending on how you choose to use it, the User guide will provide step-by-step instructions to help you complete the inversion based on the mini-LOMOS sensors. If you would like to learn more about the monitoring system, please refer to "LOMOS-mini: A coupled system quantifying transient water and heat exchanges in streambeds" by Cucchi et al. (2017). For a description of the Ginette code, please refer to Riviere et al. (2014) and Riviere et al (2019).

## Features
- Estimate water and heat fluxes by inverting temperature and pressure differential data.
- User-friendly interface for data input and output with the capacity to efficiently manage large datasets 
- The output is presented in a visually appealing and straightforward format.

## Benefits
- Understanding of the water and heat exchange processes in the hyporheic zone.
- Use the results to inform management decisions related to water resources and subsurface temperature.




## Prelude

To install gfortran and R, users can follow these steps:

Install gfortran: gfortran is a free, open-source Fortran compiler that can be used to build and execute Fortran programs. Users can install gfortran by downloading the installer from the official website (https://gcc.gnu.org/) and then following the steps given.

Install R: R is a programming language and environment for statistical computing and graphics. To install R, users can download the latest version from the official website (https://cran.r-project.org/) and follow the instructions provided to complete the installation.

#### Required_library

In order to make the R scripts work, you'll also need the packages "lubridate", "hydroGOF", "stringr" ,"stats"," stringr","ggplot2","RColorBrewer","data.table","readr", "reshape2", "akima", "cowplot," " plyr".These packages can be installed on R using the following command line:
```
install.packages(lubridate)
```

Based on your desired usage of Python or batch files and your operating system, you will need to install either Python or batch.


 


# Files preparation

## 1) Add field data.- 1) Add field data.

There must be 2 ".csv" files (one for pressure and one for temperatures) with the following format: sensor-name_point_dd_MM_YY.csv in the mini-LOMOS repertory.

#### Pay close attention to the date formatting in the titles, as this may explain why the 1D model is inoperable.Pay close attention to the date formatting in the titles, as this may explain why the 1D model is inoperable.
-------------------------------------
example (respectively for temperature and pressure):
```
t502_Point1_07_11_14.csv
p502_Point1_07_11_14.csv
```

The format of the pressure differential file is: #,dates, pressure_differential [m], temperature in stream [C]
```
1,04/11/2014 17:00:00,1.184,11.516,
etc.
```

The format of the temperature file is: #,dates, temperature depth 1 [C], temperature depth 2 [C], temperature depth 3 [C], temperature depth 4 [C]
```
1,04/11/2014 17:00:00,12.775,12.92,13.69,13.112,
etc.
```


## 2) Complete the file : inversion.COMM

An R script will use this file to read the data. 

Format the inversion.COMM command file columns as follows:
```
 $1 : point name (no space) ex: Point1
 $2 : year (two last characters YY) ex for 2014: 14 
 $3 : month (two last characters MM) ex for November: 11 
 $4 : day (two last characters DD) ex: 07 
 $5 : name temperature sensor (four first characters) ex: t502
 $6 : name pressure sensor (four first characters) ex: p502
 $7 : time step (s) ex: 900
 $8 : nb of observations ex: 3 if the 4 PT100 work, 2 if one is broken, etc. $8 = nb of working PT100-1; as the deepest working PT100 is used as a boundary condition. The upper boundary condition is the temperature sensor in the stream, in the pressure differential file.
 $9 : number of zones (clay, sand...) ex: 1 More than 2 would require more PT100 so it wouldn't make a lot of sense here.
 $10 : thickness of the upper zone (m) distance between the riverbed and the bottom of the upper zones. Type in "0" if there is just 1 zone.
 $11: maximum duration of the inversion (s) ex : 864000 (= 10 days)
```

If the duration of the field fexperiment is shorter than the prescribed duration, the shorter duration will be utilized.
------------------------------------
Example with every PT100 working and 2 zones:
```
Point1 14 11 07 t502 p504 900 3 2 0.20 864000
------------------------------------
```
Example with 2 PT100 working and 1 zone:
```
Point43 17 05 18 t520 p520 900 2 1 0 864000
```


## 3) Complete the file : inversion_PT100.COMM

Each line displays the distance between two successive working PT100s. The first line represents the space between the riverbed and the first operational PT100, the second line represents the space between the first operational PT100 and the second operational PT100, etc.
------------------------------------
example for four PT100 with 10 cm between each:
```
10
10
10
10
```
------------------------------------
example for 3 PT100 with 10 cm between each, but one is broken:
```
10
10
20
```


## 4) Completion of the inversion parameter.COMM file

Column 1 = variables of interest:

k = intrinsic permeability [m2]

n = porosity

l = solid thermal conductivity [W m−1 K−1]

r = density of solid [kg m-3]. 

The bulk volumetric heat capacity of the porous medium is calculated by the following equation :

c_mr_m = c_w r_w n + c_s r (1-n)

c_w = specific heat capacity of water [334 000 J kg−1 K−1]

r_w = density of water [1 000 kg m3]. 

c_s = solid specific heat capacity [J kg1 K1]  (imposed value) 

Column 2 = number of zones (e.g., clay, sand, etc.)

Minimum and maximum values in columns 3 and 4 indicate the range for the value test. 

Column 5 = number of tests within the previously defined range



Column 5 will reflect the number of simulations that will be launched.
------------------------------------
example for 1 zone:
```
k 1 0001D-15 0001D-11 0003
n 1 0.40 0.40 0001
l 1 1300D-03 8400D-03 0002
r 1 2600D+00 2600D+00  0001
```
------------------------------------
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
 # CHeck the input data 
The R script check_data.R can be used to verify the input data. The resulting plot will be stored in the PLOT/CHECK_DATA directory. Please be mindful to remove any outlier values before running the script.
# Userguide
Batch Userguide is available in [User_guide.md](User_guide.md)
Python Userguide is available in [User_guide_python.md](User_guide_python.md)

# OUTPUTS
The simulated results per mesh are in ginette/application/mini-LOMOS/GINETTE_SENSI/OUTPUT/

This repertory contains:


- S_flux_therm_velocity_1_t_1.dat : time, conductive heat flux, advective heat flux, total heat flux, water exhganges, temperature(1) ,temperature(2)
- Sim_heat_flux_profil_t_1.dat : time, depth,  conductive heat flux, advective heat flux, total heat flux
- Sim_temperature_maille1_1.dat: time, temperature
- Sim_temperature_maille2_1.dat:  time, temperature
- Sim_temperature_maille3_1.dat: time, temperature
- Sim_temperature_profil_t_1.dat:  time, depth , temperature
- Sim_velocity_profil_t_1.dat:  time, z, water velocity
- S_vitesse_nmaille2_hb_1.dat: time velocity top cell 1, velocity bottom cell1




"Comparaison_mailles_sim-obs.R" generates a simulated versus observed table in ginette/application/mini-LOMOS/SENSI/. It is named "Results_Stats_sim-obs".

This file contains "Sum KGE";"KGE m1";"KGE m2";"KGE m3";"RMSE m1";"RMSE m2";"RMSE m3";"MAE m1";"MAE m2";"MAE m3";"COR m1";"COR m2";"COR m3";"PBIAS m1";"PBIAS m2";"PBIAS m3" for each simulatioThe simulated results per mesh are in _ginette/application/mini-LOMOS/GINETTE_SENSI/OUTPUT/_

A comparative table of simulated versus observed is generated in _ginette/application/mini-LOMOS/SENSI/_ by "Comparaison_mailles_sim-obs.R". It is named "Results_Stats_sim-obs".


_____________________________________________________________________________________________________________________________
## PLOTS

run the script LOMOS_ginette_plot.py

The repertoire PLOT encompasses the entirety of the plot:

 The time series of  temperatures (temperature_time_series.png), the stream and aquifer water exchanges (velocity_profile.png). The total, conductive, and advective heat fluxes are in the plot (flux_timeseries.png). The interpolation of temperature data based on depth and time is depected in the plot (jolie_frise_fr.png jolie_frise_EN.png). All fluxes will be interpolated in this figures.

# License

-------------------------------------------------------------------------
This software is governed by the CeCILL v2 license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
http://www.cecill.info. 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

-------------------------------------------------------------------------
Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site http://www.cecill.info.

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.


