#!/bin/bash
# Quick bash script that puts together all relevant measurements

DIR="/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/EMD_Data"

for day in $DIR/* 
do
    mkdir "$day"/0_Compressed --p
    # TODO : Separate in different Regions for comparison
    for is_kind in "Solo_N" "Solo_T" "Solo_V_R" "Solo_Mf" ; do
        for PERIOD in "3 - 20"; do
            # Do for all wavelengths separately
            for WV in 94 193 211; do
                mkdir "$day"/0_Compressed/"$PERIOD"/"$WV" --p
                cp "$day"/*/$is_kind/*/"$PERIOD"/*"$WV"*.png "$day"/0_Compressed/"$PERIOD"/"$WV"/
            done
        done
    done
done    

echo "DONE in $DIR"