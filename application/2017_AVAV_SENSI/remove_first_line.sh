awk -f remove_first_line.awk $1
echo $2
mv awk.out tmp$2
#rm $1
