This User_guide will guide you step by step to help you complete the inversion based on the mini-LOMOS sensors, depending on how you chose to use it. 
Please refer to "LOMOS-mini: A coupled system quantifying transient water and heat exchanges in streambeds" from Cucchi et al. (2017) if you want to know more information about the monitoring system and to  "Experimental and numerical assessment of transient stream-aquifer exchange during disconnection" from Riviere et al. (2014) and ."Pore water pressure evolution below a freezing front under saturated conditions: Large-scale laboratory experiment and numerical investigation." Rivière et al. (2019) for the numerical model.
_____________________________________________________________________________________________________________________________

### Prelude :

Before making anything, please make sure you've installed the "dos2unix" UNIX package. You can do it from the terminal:
> apt install dos2unix

Also, make sure that all the files of the "mini-LOMOS" file are executable:
> chmod 775 *______________________________________________________________________________________________________________________

There is nothing for Unix users to do.
### windows-users
To install sed and awk on a Windows machine, you can follow these steps:

- 1) Download the sed and awk executables: You can find pre-built sed and awk executables for Windows on the GnuWin32 website (https://sourceforge.net/projects/gnuwin32/files/sed/). Download the executables and save them to a directory on your machine.

- 2) Add the directory to your PATH environment variable: To use sed and awk from the command line, you will need to add the directory where you saved the executables to your PATH environment variable. To do this, follow these steps:
  - Press the Windows key and type "Environment Variables."
  - Click on the "Edit the system environment variables" button.
  - In the "System Properties" window, click on the "Environment Variables" button.
  - In the "Environment Variables" window, scroll down to the "System Variables" section and find the "Path" variable.
  - Click on the "Edit" button.
  - In the "Edit environment variable" window, click on the "New" button and add the path to the directory where you saved the sed and awk executables.
  - Click "OK" to close all windows.

- 3) Test the installation: To test the installation, open a command prompt and type sed --version or awk --version. If the executables are installed correctly, you should see the version number displayed.

_
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
