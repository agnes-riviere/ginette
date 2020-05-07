#!/bin/bash
awk -f remove_first_line.awk $1
mv awk.out tmp$2
#rm $1
