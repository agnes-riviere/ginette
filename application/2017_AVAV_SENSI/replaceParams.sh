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


# replace intrinsic permeability field [k] in E_zone_parameter.dat by first argument
sed -i -e 's/\[k2\]/'$3'/' E_zone_parameter.dat

# replace porosity field [n] in E_zone_parameter.dat by second argument
sed -i -e 's/\[n2\]/'$4'/' E_zone_parameter.dat


# replace intrinsic permeability field [k] in E_zone_parameter.dat by first argument
sed -i -e 's/\[k3\]/'${5}'/' E_zone_parameter.dat

# replace porosity field [n] in E_zone_parameter.dat by second argument
sed -i -e 's/\[n3\]/'${6}'/' E_zone_parameter.dat
# replace intrinsic permeability field [k] in E_zone_parameter.dat by first argument
sed -i -e 's/\[k4\]/'${7}'/' E_zone_parameter.dat

# replace porosity field [n] in E_zone_parameter.dat by second argument
sed -i -e 's/\[n4\]/'${8}'/' E_zone_parameter.dat

# replace intrinsic permeability field [k] in E_zone_parameter.dat by first argument
sed -i -e 's/\[k5\]/'${9}'/' E_zone_parameter.dat

# replace porosity field [n] in E_zone_parameter.dat by second argument
sed -i -e 's/\[n5\]/'${10}'/' E_zone_parameter.dat

# replace intrinsic permeability field [k] in E_zone_parameter.dat by first argument
sed -i -e 's/\[k6\]/'${11}'/' E_zone_parameter.dat

# replace porosity field [n] in E_zone_parameter.dat by second argument
sed -i -e 's/\[n6\]/'${12}'/' E_zone_parameter.dat

# replace intrinsic permeability field [k] in E_zone_parameter.dat by first argument
sed -i -e 's/\[k7\]/'${13}'/' E_zone_parameter.dat

# replace porosity field [n] in E_zone_parameter.dat by second argument
sed -i -e 's/\[n7\]/'${14}'/' E_zone_parameter.dat


# replace intrinsic permeability field [k] in E_zone_parameter.dat by first argument
sed -i -e 's/\[k8\]/'${15}'/' E_zone_parameter.dat

# replace porosity field [n] in E_zone_parameter.dat by second argument
sed -i -e 's/\[n8\]/'${16}'/' E_zone_parameter.dat


