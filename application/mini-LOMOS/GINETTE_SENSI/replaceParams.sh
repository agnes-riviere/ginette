#!/bin/bash
# taking as argument the numerical value for parameters in the loop
# taking as argument the numerical value for parameters in the loop
# There are 3 arguments, in order :
#  akx and akz  tagged by [k]
# thermal conductivity - tagged by [l]
#  porosity -  tagged by [n]
#- solid thermal density  by [r]
# replacing tags eg.[n] by numerical values in ginette input files




# replace intrinsic permeability field [k] in E_zone_parameter.dat by first argument
sed -i -e 's/\[k1\]/'$1'/' E_zone_parameter.dat

# replace porosity field [n] in E_zone_parameter.dat by second argument
sed -i -e 's/\[n1\]/'$2'/' E_zone_parameter.dat

# replace thermal conductivity [n] in E_zone_parameter.dat by third argument
sed -i -e 's/\[l1\]/'$3'/' E_zone_parameter.dat

# replace solid thermal density[rho_s] in E_zone_parameter.dat by fourth argument
sed -i -e 's/\[r1\]/'$4'/' E_zone_parameter.dat

