#!/bin/bash


#### print title

CORAL='\033[38;5;208m'
GOLD='\033[38;5;220m'
CYAN='\033[36m'
RED='\033[31m'
GREEN='\033[38;5;40m'

RESET='\033[0m'
BOLD='\033[1m'


echo -e "${GOLD}${BOLD}"
echo "  ███████╗██████╗  █████╗  ██████╗████████╗ █████╗ ██╗     ██████╗ ███████╗███████╗██████╗ "
echo "  ██╔════╝██╔══██╗██╔══██╗██╔════╝╚══██╔══╝██╔══██╗██║     ██╔══██╗██╔════╝██╔════╝██╔══██╗"
echo "  █████╗  ██████╔╝███████║██║        ██║   ███████║██║     ██║  ██║█████╗  █████╗  ██████╔╝"
echo "  ██╔══╝  ██╔══██╗██╔══██║██║        ██║   ██╔══██║██║     ██║  ██║██╔══╝  ██╔══╝  ██╔══██╗"
echo "  ██║     ██║  ██║██║  ██║╚██████╗   ██║   ██║  ██║███████╗██████╔╝███████╗███████╗██║  ██║"
echo "  ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝   ╚═╝   ╚═╝  ╚═╝╚══════╝╚═════╝ ╚══════╝╚══════╝╚═╝  ╚═╝"
echo -e "${RESET}"

# Sous-titre et informations du projet
echo -e "  ${CYAN}${BOLD}══ L1 PIPELINE RUNNER ═════════════════════════════════════════════════════════${RESET}"
echo -e "  ${BOLD}Version    :${RESET} 1.0"
echo -e "  ${BOLD}Programmeur:${RESET} ${CORAL}Lucas Villard / LuVil86${RESET}"
echo -e "  ${BOLD}Date       :${RESET} $(date +'%d/%m/%Y à %H:%M')"
echo -e "  ${CYAN}${BOLD}═══════════════════════════════════════════════════════════════════════════════════${RESET}"
echo "" 


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
if [[ "${currentEnvir}" != "pipeline_fractalDeer_prd" ]]
then
echo -e "${RED} ERROR ! current environment is not 'pipeline_fractalDeer_prd' : aborting execution"
echo -e "${RED} - in case you haven't created the environment, run the 'environment_fractalDeer.yml' file in conda "
echo -e "${RED}   with the command 'conda env create -f ./Workflow_fractalDeer/PIPELINE_FRACTALDEER/environment_fractalDeer_PRD.yml' "
echo -e "${RED}   and activate it using the commmand 'conda activate pipeline_fractalDeer_prd'${RESET}"

exit 2
fi



inPath=$(realpath $1)
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


echo -e "${GREEN}${BOLD}input csv file --> ${inPath}${RESET}" 

echo -e "${GREEN}${BOLD}outFolder -->  ${outFolder}${RESET}"

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


if [[ $? -eq 1 ]]
then
echo " ************** Script S1 returned an error : aborting pipeline ************************"
exit 2   
fi

### -2. ** compute fractal statistics ** #####
python ./S2_compute_Fractal_observed_and_simulated_trajectory.py "${inPath}" "${outFolder}"


if [[ $? -eq 1 ]]
then
echo " ************** Script S2 returned an error : aborting pipeline ************************"
exit 2   
fi

### -3. ** obtain Zscores ** ####
/usr/bin/Rscript ./S3_get_Zscores_CRW.R "${outFolder}"

if [[ $? -eq 1 ]]
then
echo " ************** Script S3 returned an error : aborting pipeline ************************"
exit 2   
fi
### -3.1 : if list of candidates is empty, abort pipeline #######
stepList=$(find "${outFolder}" -regex ".*_candidate_.*")

dimList=$(echo "${stepList}" | grep -c '^')

# 3. Si aucun fichier n'est trouvé (dimList -eq 0)
if [[ $dimList -eq 0 ]]; then
    echo "the step Size list is empty : aborting the next steps of the pipeline"
    exit 0
fi


echo "Candidate stepSizes have been found ! reading : " "${stepList}"


### -4. ** generate behaviour classification from candidate list ** ###
# Attention : si ton script Python attend un seul fichier mais que 'find' en a trouvé plusieurs,
# il faudra adapter la boucle. Ici on passe la liste brute.
python ./S4_get_behaviour_vector_from_list.py "${inPath}" "${stepList}"
#~ ### -5. ** compute Chi-square on a user-specified path metric ** ####
#~ /usr/bin/Rscript ./S5_choose_best_stepSize.R "${outFolder}" "${stepList}" "ratio_endNSD"
#~ /usr/bin/Rscript ./S5_choose_best_stepSize.R "${outFolder}" "${stepList}" "ratio_meanNSD"
#~ /usr/bin/Rscript ./S5_choose_best_stepSize.R "${outFolder}" "${stepList}" "ratio_cumulativeNSD"

if [[ $? -eq 1 ]]
then
echo " ************** Script S4 returned an error : aborting pipeline ************************"
exit 2   
fi
/usr/bin/Rscript ./S5_choose_best_stepSize_V3.R "${outFolder}" "${stepList}"

if [[ $? -eq 1 ]]
then
echo " ************** Script S5 returned an error : aborting pipeline ************************"
exit 2   
fi

end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"


# Titre géant en ASCII Art : SUCCESS
echo -e "${GREEN}${BOLD}"
echo "  ███████╗██╗   ██╗ ██████╗ ██████╗███████╗███████╗███████╗"
echo "  ██╔════╝██║   ██║██╔════╝██╔════╝██╔════╝██╔════╝██╔════╝"
echo "  ███████╗██║   ██║██║     ██║     █████╗  ███████╗███████╗"
echo "  ╚════██║██║   ██║██║     ██║     ██╔══╝  ╚════██║╚════██║"
echo "  ███████║╚██████╔╝╚██████╗╚██████╗███████╗███████║███████║"
echo "  ╚══════╝ ╚═════╝  ╚═════╝ ╚═════╝╚══════╝╚══════╝╚══════╝"
echo -e "${RESET}"

# Sous-titre et résumé de clôture
echo -e "  ${CYAN}${BOLD}══ PROCESS TERMINATED ═════════════════════════════════════${RESET}"
echo -e "  ${BOLD}Statut :${RESET} ${GREEN}${BOLD}PIPELINE ENDED SUCCESSFULLY${RESET}"
echo -e "  ${BOLD}Heure  :${RESET} $(date +'%H:%M:%S')"
echo -e "Processing finished in $elapsed seconds"
echo -e "  ${CYAN}${BOLD}═══════════════════════════════════════════════════════════${RESET}"
echo "" # Ligne vide finale
