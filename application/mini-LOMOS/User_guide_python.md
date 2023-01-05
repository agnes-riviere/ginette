This User_guide will guide you step by step to help you complete the inversion based on the mini-LOMOS sensors, depending on how you chose to use it. 
Please refer to "LOMOS-mini: A coupled system quantifying transient water and heat exchanges in streambeds" from Cucchi et al. (2017) if you want to know more information about the monitoring system and to  "Experimental and numerical assessment of transient stream-aquifer exchange during disconnection" from Riviere et al. (2014) and ."Pore water pressure evolution below a freezing front under saturated conditions: Large-scale laboratory experiment and numerical investigation." Rivière et al. (2019) for the numerical model.
_____________________________________________________________________________________________________________________________

### Prelude :
To Install Python: Python is a popular, high-level programming language that can be used for a wide range of tasks. Users can install Python by downloading the most recent version from the official website (https://www.python.org/) and then following the steps given.

#### required_modules

os, sys, numpy, pathlib, scipy, matplotlib, subprocess,rpy2, pandas

__________________________________________________________________________________________________________________
### 1) Launch the screening of multiple experiments of temperature and pressure in the HZ with the awk script inversion.awk from the terminal:
```
> awk -f inversion.awk inversion.COMM
```
### 2) To run a simulation using nohup, you can use the following command:
```
> nohup awk -f inversion.awk inversion.COMM &
```
This will run the simulation script in the background, even if you close the terminal window. The nohup command stands for "no hangup," and it prevents the process from being terminated when the terminal is closed. The & symbol at the end of the command causes the process to run in the background.
_________________________________________________________________________________________________________________________
### 3) OUTPUT
The simulated results per mesh are in _ginette/application/mini-LOMOS/GINETTE_SENSI/OUTPUT/_

A comparative table of simulated versus observed is generated in _ginette/application/mini-LOMOS/SENSI/_ by "Comparaison_mailles_sim-obs.R". It is named "Results_Stats_sim-obs".

**Don't forget to copy/paste it in a specific repository if you want to keep track of these files, as launching another simulation will overwrite them.**
_____________________________________________________________________________________________________________________________
### 7) PLOTS

**Results:**

When you achieved satisfying calibrations, you can plot the simulated vs observed temperature time-series thanks to the temperature_time_series.R script in ginette/application/mini-LOMOS/PLOT/
