#!/bin/bash
# Quick bash script that puts together all relevant measurements

DIR="/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/EMD_Data"
echo $DIR
mkdir $DIR --p
for day in "AIA_D28_H4:5" "AIA_D28_H5:6" "AIA_D28_H6:7"
do
    mkdir $DIR/$day/Compressed --p
    for rs_kind in "Lcurve_16" "Lcurve_17" "Lcurve_18" "Lcurve_21" "Lcurve_22" "Lcurve_23"
        do
        for is_kind in "Solo_N"
            do 
            cp $DIR/$day/$rs_kind/$is_kind/*/*/*$rs_kind*.png $DIR/$day/Compressed/
            done
        done
done    
