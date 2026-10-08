#!/bin/bash

cd ~/Workflow_fractalDeer/PIPELINE_FRACTALDEER

currentEnvir=$(echo $CONDA_PREFIX | awk -F "/" '{print $NF}')
echo "current Environment : " ${currentEnvir}
if [[ "${currentEnvir}" != "pipeline_fractalDeer" ]]
then
echo "current environment is not 'pipeline_fractalDeer' : aborting execution"
echo " - in case you haven't created the environment, run the 'environment_fractalDeer.yml' file in conda "
echo "   with the command 'conda env create -f /home/luvil/Workflow_fractalDeer/PIPELINE_FRACTALDEER/environment_fractalDeer.yml' "
echo "   and activate it using the commmand 'conda activate pipeline_fractalDeer'"

exit 2
fi

#### please put the full path (not the ~ alias) + the folders should exists 
inFolder="/home/luvil/IE-OFEV/cerf_grisons_3"
outFolder="/home/luvil/IE-OFEV/results_grisons_3"
deerList="/home/luvil/IE-OFEV/cerf_grisons_3/RED_DEER_GRISONS_3_70_percent.csv"


while IFS=, read -r prenom deerYear nbFixes	startDT endDT Duration_days	is_full_year
do


inFile=${inFolder}/${prenom}/${prenom}_${deerYear}.csv

echo " %%  input file -->> ${inFile}" 


if eval ! test -e "$inFile"
then
echo "LAUNCHER- The input file trajectory does not exists :: aborting launcher"
exit 2 
fi


eval /home/luvil/Workflow_fractalDeer/PIPELINE_FRACTALDEER/L1_launch_fractalDeer.sh ${inFile} ${outFolder}
if [[ $? -eq 2 ]]
then
echo "LAUNCHER- One of the pipeline in the list failed :: aborting launcher"
exit 2
fi
done < <(tail -n+2 $deerList)



