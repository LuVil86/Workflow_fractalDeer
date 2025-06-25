#!/bin/bash

cd ~/Workflow_fractalDeer/PIPELINE_FRACTALDEER
outFolder="~/results_grisons/"
while IFS=, read -r prenom deerYear nbFixes	startDT endDT Duration_days	is_full_year
do


inFile=~/IE-OFEV/grisons_data/red_deer_grisons_parsed/${prenom}/${prenom}_${deerYear}.csv
eval /home/luvil/Workflow_fractalDeer/PIPELINE_FRACTALDEER/L1_launch_fractalDeer.sh ${inFile} ${outFolder}
done < <(tail -n+2 ~/IE-OFEV/grisons_data/subset_red_deer_70percent_days.csv)
