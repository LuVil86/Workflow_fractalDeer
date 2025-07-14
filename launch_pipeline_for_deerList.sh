#!/bin/bash

cd ~/Workflow_fractalDeer/PIPELINE_FRACTALDEER

outFolder="~/IE-OFEV/results_jura"
while IFS=, read -r prenom deerYear nbFixes	startDT endDT Duration_days	is_full_year
do


inFile=~/IE-OFEV/cerf_jura/${prenom}/${prenom}_${deerYear}.csv
eval /home/luvil/Workflow_fractalDeer/PIPELINE_FRACTALDEER/L1_launch_fractalDeer.sh ${inFile} ${outFolder}
if [[ $? -eq 2 ]]
then
echo "LAUNCHER- One of the pipeline in the list failed :: aborting launcher"
exit 2
fi
done < <(tail -n+2 /home/luvil/IE-OFEV/cerf_jura/RED_DEER_JURA_fixes_counts_subset_70_percent.csv)
