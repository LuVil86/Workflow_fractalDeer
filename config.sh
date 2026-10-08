#!/bin/bash
#conda env update --file ./environment_fractalDeer_PRD.yml --prune
conda env create --file ./environment_fractalDeer_PRD.yml 
~/anaconda3/envs/pipeline_fractalDeer_prd/bin/python -m ipykernel install --user --name pipeline_fractalDeer_prd --display-name "Conda (pipeline_fractalDeer)"

echo "Conda environement successfully installed"
