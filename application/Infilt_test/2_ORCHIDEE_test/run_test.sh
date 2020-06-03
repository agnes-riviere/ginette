#!/bin/bash

DIRTEST=/home/dkilic/Modèles/ginette3/application/Warrick/Special_cases/2_ORCHIDEE_test

topBname=E_debit_haut_t.dat_veg
tobBfinal=E_debit_haut_t.dat

cd $DIRTEST

for i in {5..15}
do
    cp E_debit_haut_t.dat_veg$i $tobBfinal
    ./ginette
    wait
    mv S_* out_veg_$i
done

