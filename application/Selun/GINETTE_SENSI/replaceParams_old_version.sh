#!/bin/bash
# taking as argument the numerical value for parameters in the loop
# There are 3 arguments, in order :
# $1 - intrinsic permeability - akx and akz in E_parametre.dat - tagged by [k]
# $2 - solid thermal conductivity - alandam in E_p_therm.dat - tagged by [lambda_s]
# $3 - porosity - omp in E_parametre.dat - tagged by [n]
# $4 - solid thermal capacity - cpm in E_p_therm.dat - tagged by [c_s]
# $5 - solid thermal density - rhosi in E_p_therm.dat - tagged by [rho_s]
# replacing tags eg.[n] by numerical values in ginette input files

# replace intrinsic permeability field [k] in E_parametre.dat by first argument
sed -i -e 's/\[k\]/'$1'/' E_parametre.dat

# replace porosity field [n] in E_parametre.dat by second argument
sed -i -e 's/\[n\]/'$2'/' E_parametre.dat

# replace porosity field [n] in E_parametre.dat by third argument
sed -i -e 's/\[sto\]/'$3'/' E_parametre.dat

# replace solid thermal capacity [c_s] in E_p_therm.dat by third argument
sed -i -e 's/\[lambda_s\]/'$4'/' E_p_therm.dat

# replace solid thermal density [rho_s] in E_p_therm.dat by fourth  argument
sed -i -e 's/\[rho_s\]/'$5'/' E_p_therm.dat

