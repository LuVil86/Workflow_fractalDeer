#!/bin/bash


if [[ -z "$1" ]]
then
echo "you must specify a input CSV of GPS data!!"
exit 2
fi


if [[ -z "$2" ]]
then
echo "you must specify an output folder location !!"
exit 2
fi

#### environement check #####

cd "$(dirname "$0")"

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



inPath=$( realpath $1)
inCSVfile=$(basename $1 .csv)
#LOCALIZE=$(echo ${inPath%/*} | sed 's/[\\ '$'\t'']/\\&/')
deerName=$(echo $inCSVfile | cut -d "_" -f 1)
deerYear=$(echo $inCSVfile| cut -d "_" -f 2,3)

## version if outFolder is directly related to inFile
#echo ${LOCALIZE}
#LOCALIZE=${inPath%/*}
#outFolder=${LOCALIZE}/OUT_${deerName}_${deerYear}

### version where we specified the path where we want the outfolder to be into
outFolder=$(realpath $2)/OUT_${deerName}_${deerYear}

echo "-- L1 PIPELINE FRACTALDEER V 1.0 --"

echo "input csv file --> " ${inPath}

echo "outFolder --> " ${outFolder}
if [[ ! -d "${outFolder}" ]] 
then
mkdir "${outFolder}" 

fi

start_time="$(date -u +%s)"

##### -0. plot migrateR results ####
/usr/bin/Rscript ./S0_plotMigrateR.R "${inPath}" "${outFolder}" TRUE
if [[ $? -eq 1 ]]
then
echo " ************** Script S0 returned an error : aborting pipeline ************************"
exit 2   
fi
### -1. ** generate CRW ** ####
/usr/bin/Rscript ./S1_generate_CRW.R "${inPath}" "${outFolder}"

### -2. ** compute fractal statistics ** #####
python ./S2_compute_Fractal_observed_and_simulated_trajectory.py "${inPath}" "${outFolder}"

### -3. ** obtain Zscores ** ####
/usr/bin/Rscript ./S3_get_Zscores_CRW.R "${outFolder}"

### -3.1 : if list of candidates is empty, abort pipeline #######
stepList=$(find "${outFolder}" -regex ".*_candidate_.*")
echo "${stepList}"
dimList=$(wc -l "${stepList}" | cut -d " " -f 1)
if [[ $dimList -le 1 ]];
then
echo "the step Size list is empty : aborting the next steps of the pipeline"
exit 0
fi
### -4. ** generate behaviour classification from candidate list ** ###
python ./S4_get_behaviour_vector_from_list.py "${inPath}" "${stepList}"
#~ ### -5. ** compute Chi-square on a user-specified path metric ** ####
#~ /usr/bin/Rscript ./S5_choose_best_stepSize.R "${outFolder}" "${stepList}" "ratio_endNSD"
#~ /usr/bin/Rscript ./S5_choose_best_stepSize.R "${outFolder}" "${stepList}" "ratio_meanNSD"
#~ /usr/bin/Rscript ./S5_choose_best_stepSize.R "${outFolder}" "${stepList}" "ratio_cumulativeNSD"

/usr/bin/Rscript ./S5_choose_best_stepSize_V3.R "${outFolder}" "${stepList}"

echo "--- PIPELINE ENDED SUCCESSFULLY -- "
end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"

echo "Processing finished in $elapsed seconds"

