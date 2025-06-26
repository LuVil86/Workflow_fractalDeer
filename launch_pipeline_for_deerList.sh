#!/bin/bash

cd ~/Workflow_fractalDeer/PIPELINE_FRACTALDEER

outFolder="/media/luvil/T7/IE-OFEV/results_Jura"
while IFS=, read -r prenom deerYear nbFixes	startDT endDT Duration_days	is_full_year
do


inFile=/media/luvil/T7/IE-OFEV/cerf_jura/${prenom}/${prenom}_${deerYear}.csv
eval /home/luvil/Workflow_fractalDeer/PIPELINE_FRACTALDEER/L1_launch_fractalDeer.sh ${inFile} ${outFolder}
done < <(tail -n+2 /media/luvil/T7/IE-OFEV/cerf_jura/RED_DEER_JURA_fixes_counts_subset_70_percent.csv)
