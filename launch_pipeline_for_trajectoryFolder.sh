#!/bin/bash

cd ~/Workflow_fractalDeer/PIPELINE_FRACTALDEER

inFolder="/home/luvil/IE-OFEV/new_cerfs_valais_2025/trajectories_denoised/"
outFolder="~/IE-OFEV/results_Valais"

for inFile in $(ls ${inFolder})
do
	echo ${inFile}
	eval /home/luvil/Workflow_fractalDeer/PIPELINE_FRACTALDEER/L1_launch_fractalDeer.sh ${inFolder}${inFile} ${outFolder}
	if [[ $? -eq 2 ]]
	then
		echo "LAUNCHER- One of the pipeline in the list failed :: aborting launcher"
	exit 2
	fi
done 
