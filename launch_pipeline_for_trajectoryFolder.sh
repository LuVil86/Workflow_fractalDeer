#!/bin/bash

if [[ -z $1 ]]
then
    echo "you must specify an input folder where the trajectories are located !!"
    exit 2
fi

if [[ -z $2 ]]
then
    echo "you must specify an output folder location !!"
    exit 2
fi

inFolder=$1
outFolder=$2

if [[ ! -d $inFolder ]]
then
    echo "the input folder does not exist : aborting script"
    exit 2
fi

if [[ -z $(shopt -s nullglob; echo "$inFolder"/*) ]]
then
    echo "the input folder is empty : aborting script"
    exit 2
fi

if [[ ! -d $outFolder ]]
then
    echo "the output folder specified does not exist : creating it.."
    mkdir -p "$outFolder"
fi

for inFilePath in "$inFolder"/*
do
    inFile=$(basename "$inFilePath")
    
    echo "${inFile}"
    
    ~/Workflow_fractalDeer/PIPELINE_FRACTALDEER/L1_launch_fractalDeer.sh "${inFolder}/${inFile}" "${outFolder}"
    
    if [[ $? -eq 2 ]]
    then
        echo "LAUNCHER- One of the pipeline in the list failed :: aborting launcher"
        exit 2
    fi
done
