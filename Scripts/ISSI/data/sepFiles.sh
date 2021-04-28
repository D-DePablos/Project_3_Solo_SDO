 #!/bin/bash
# Quick bash script that puts together all relevant measurements

DIR="/home/diegodp/Documents/PhD/Paper_3/SolO_SDO_EUI/unsafe/ISSI"

for day in $DIR/* 
do
    mkdir "$day"/Compressed --p
    # TODO : Separate in different Regions for comparison
    for is_kind in "171_plume" "???" 
    do
        for PERIOD in "5 - 60"; do
            mkdir "$day"/Compressed/"$PERIOD"/ --p
            cp "$day"/*/$is_kind/*/"$PERIOD"/*.png "$day"/Compressed/"$PERIOD"/
        done
    done
done    

echo "DONE in $DIR"
